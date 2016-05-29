#include "miner.h"
#include "algo-gate-api.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "algo/groestl/sph_groestl.h"
#include "algo/blake/sph_blake.h"
#include "algo/bmw/sph_bmw.h"
#include "algo/jh/sph_jh.h"
#include "algo/keccak/sph_keccak.h"
#include "algo/skein/sph_skein.h"
#include "algo/shavite/sph_shavite.h"
#include "algo/luffa/sph_luffa.h"
#include "algo/cubehash/sph_cubehash.h"
#include "algo/simd/sph_simd.h"
#include "algo/echo/sph_echo.h"
#include "algo/hamsi/sph_hamsi.h"
#include "algo/fugue/sph_fugue.h"

#include "algo/luffa/sse2/luffa_for_sse2.h" 
#include "algo/cubehash/sse2/cubehash_sse2.h"
#include "algo/simd/sse2/nist.h"
#include "algo/blake/sse2/blake.c"   
#include "algo/bmw/sse2/bmw.c"
#include "algo/keccak/sse2/keccak.c"
#include "algo/skein/sse2/skein.c"
#include "algo/jh/sse2/jh_sse2_opt64.h"

#ifdef NO_AES_NI
  #include "algo/groestl/sse2/grso.h"
  #include "algo/groestl/sse2/grso-macro.c"
#else
  #include "algo/groestl/aes_ni/hash-groestl.h"
  #include "algo/echo/aes_ni/hash_api.h"
#endif

typedef struct {
#ifdef NO_AES_NI
        sph_groestl512_context   groestl;
        sph_echo512_context      echo;
#else
        hashState_groestl       groestl;
        hashState_echo          echo;
#endif
        hashState_luffa         luffa;
        cubehashParam           cubehash;
        sph_shavite512_context  shavite;
        hashState_sd            simd;
        sph_hamsi512_context    hamsi;
        sph_fugue512_context    fugue;
} x13_ctx_holder;

x13_ctx_holder x13_ctx;

void init_x13_ctx()
{
#ifdef NO_AES_NI
        sph_groestl512_init(&x13_ctx.groestl);
        sph_echo512_init(&x13_ctx.echo);
#else
        init_echo( &x13_ctx.echo, 512 );
        init_groestl (&x13_ctx.groestl );
#endif
        init_luffa( &x13_ctx.luffa, 512 );
        cubehashInit( &x13_ctx.cubehash, 512, 16, 32 );
        sph_shavite512_init( &x13_ctx.shavite );
        init_sd( &x13_ctx.simd, 512 );
        sph_hamsi512_init( &x13_ctx.hamsi );
        sph_fugue512_init( &x13_ctx.fugue );
};

static void x13hash(void *output, const void *input)
{
	unsigned char hash[128]; // uint32_t hashA[16], hashB[16];
	#define hashB hash+64
      
        x13_ctx_holder ctx;
        memcpy( &ctx, &x13_ctx, sizeof(x13_ctx) );
#ifdef NO_AES_NI
        grsoState sts_grs;
#endif

        // X11 algos

        unsigned char hashbuf[128];
        size_t hashptr;
        sph_u64 hashctA;
        sph_u64 hashctB;

        //---blake1---

        DECL_BLK;
        BLK_I;
        BLK_W;
        BLK_C;

        //---bmw2---

        DECL_BMW;
        BMW_I;
        BMW_U;

        #define M(x)    sph_dec64le_aligned(data + 8 * (x))
        #define H(x)    (h[x])
        #define dH(x)   (dh[x])

        BMW_C;

        #undef M
        #undef H
        #undef dH
        
        //---groetl----

#ifdef NO_AES_NI
// use GRS if possible
          GRS_I;
          GRS_U;
          GRS_C;
//        sph_groestl512 (&ctx.groestl, hash, 64);
//        sph_groestl512_close(&ctx.groestl, hash);
#else
        update_groestl( &ctx.groestl, (char*)hash,512);
        final_groestl( &ctx.groestl, (char*)hash);
#endif

        //---skein4---

        DECL_SKN;
        SKN_I;
        SKN_U;
        SKN_C;

        //---jh5------

        DECL_JH;
        JH_H;

        //---keccak6---

        DECL_KEC;
        KEC_I;
        KEC_U;
        KEC_C;

        //--- luffa7
        update_luffa( &ctx.luffa, (const BitSequence*)hash,512);
        final_luffa( &ctx.luffa, (BitSequence*)hashB);

        // 8 Cube
        cubehashUpdate( &ctx.cubehash, (const byte*) hashB,64);
        cubehashDigest( &ctx.cubehash, (byte*)hash);

        // 9 Shavite
        sph_shavite512( &ctx.shavite, hash, 64);
        sph_shavite512_close( &ctx.shavite, hashB);

        // 10 Simd
        update_sd( &ctx.simd, (const BitSequence *)hashB,512);
        final_sd( &ctx.simd, (BitSequence *)hash);

        //11---echo---

#ifdef NO_AES_NI
        sph_echo512(&ctx.echo, hash, 64);
        sph_echo512_close(&ctx.echo, hashB);
#else
        update_echo ( &ctx.echo, (const BitSequence *) hash, 512);
        final_echo( &ctx.echo, (BitSequence *) hashB);
#endif

        // X13 algos
        // 12 Hamsi
	sph_hamsi512(&ctx.hamsi, hashB, 64);
	sph_hamsi512_close(&ctx.hamsi, hash);

        // 13 Fugue
	sph_fugue512(&ctx.fugue, hash, 64);
	sph_fugue512_close(&ctx.fugue, hashB);

#if defined(_WIN32)	|| defined(_WIN64) || defined(__TOS_WIN__) || defined(__WINDOWS__)
    	__asm__ volatile("emms");
#else
		asm volatile ("emms");
#endif		
	memcpy(output, hashB, 32);
}

void x13hash_alt(void *output, const void *input)
{
        unsigned char hash[128]; // uint32_t hashA[16], hashB[16];
        #define hashB hash+64

        sph_blake512_context     ctx_blake;
        sph_bmw512_context       ctx_bmw;
        sph_groestl512_context   ctx_groestl;
        sph_jh512_context        ctx_jh;
        sph_keccak512_context    ctx_keccak;
        sph_skein512_context     ctx_skein;
        sph_luffa512_context     ctx_luffa;
        sph_cubehash512_context  ctx_cubehash;
        sph_shavite512_context   ctx_shavite;
        sph_simd512_context      ctx_simd;
        sph_echo512_context      ctx_echo;
        sph_hamsi512_context     ctx_hamsi;
        sph_fugue512_context     ctx_fugue;

        sph_blake512_init(&ctx_blake);
        sph_blake512(&ctx_blake, input, 80);
        sph_blake512_close(&ctx_blake, hash);

        sph_bmw512_init(&ctx_bmw);
        sph_bmw512(&ctx_bmw, hash, 64);
        sph_bmw512_close(&ctx_bmw, hashB);

        sph_groestl512_init(&ctx_groestl);
        sph_groestl512(&ctx_groestl, hashB, 64);
        sph_groestl512_close(&ctx_groestl, hash);

        sph_skein512_init(&ctx_skein);
        sph_skein512(&ctx_skein, hash, 64);
        sph_skein512_close(&ctx_skein, hashB);

        sph_jh512_init(&ctx_jh);
        sph_jh512(&ctx_jh, hashB, 64);
        sph_jh512_close(&ctx_jh, hash);

        sph_keccak512_init(&ctx_keccak);
        sph_keccak512(&ctx_keccak, hash, 64);
        sph_keccak512_close(&ctx_keccak, hashB);

        sph_luffa512_init(&ctx_luffa);
        sph_luffa512(&ctx_luffa, hashB, 64);
        sph_luffa512_close(&ctx_luffa, hash);

        sph_cubehash512_init(&ctx_cubehash);
        sph_cubehash512(&ctx_cubehash, hash, 64);
        sph_cubehash512_close(&ctx_cubehash, hashB);

        sph_shavite512_init(&ctx_shavite);
        sph_shavite512(&ctx_shavite, hashB, 64);
        sph_shavite512_close(&ctx_shavite, hash);
        sph_simd512_init(&ctx_simd);
        sph_simd512(&ctx_simd, hash, 64);
        sph_simd512_close(&ctx_simd, hashB);

        sph_echo512_init(&ctx_echo);
        sph_echo512(&ctx_echo, hashB, 64);
        sph_echo512_close(&ctx_echo, hash);
        sph_hamsi512_init(&ctx_hamsi);
        sph_hamsi512(&ctx_hamsi, hash, 64);
        sph_hamsi512_close(&ctx_hamsi, hashB);

        sph_fugue512_init(&ctx_fugue);
        sph_fugue512(&ctx_fugue, hashB, 64);
        sph_fugue512_close(&ctx_fugue, hash);

        memcpy(output, hash, 32);
}


int scanhash_x13(int thr_id, struct work *work, uint32_t max_nonce,
                               uint64_t *hashes_done)
{
        uint32_t *pdata = work->data;
        uint32_t *ptarget = work->target;

	uint32_t n = pdata[19] - 1;
	const uint32_t first_nonce = pdata[19];
	const uint32_t Htarg = ptarget[7];

	uint32_t _ALIGN(32) hash64[8];
	uint32_t endiandata[32];

	uint64_t htmax[] = {
		0,
		0xF,
		0xFF,
		0xFFF,
		0xFFFF,
		0x10000000
	};
	uint32_t masks[] = {
		0xFFFFFFFF,
		0xFFFFFFF0,
		0xFFFFFF00,
		0xFFFFF000,
		0xFFFF0000,
		0
	};

//        init_13_contexts();

	// we need bigendian data...


                be32enc( &endiandata[0], ((uint32_t*)pdata)[0] );
                be32enc( &endiandata[1], ((uint32_t*)pdata)[1] );
                be32enc( &endiandata[2], ((uint32_t*)pdata)[2] );
                be32enc( &endiandata[3], ((uint32_t*)pdata)[3] );
                be32enc( &endiandata[4], ((uint32_t*)pdata)[4] );
                be32enc( &endiandata[5], ((uint32_t*)pdata)[5] );
                be32enc( &endiandata[6], ((uint32_t*)pdata)[6] );
                be32enc( &endiandata[7], ((uint32_t*)pdata)[7] );
                be32enc( &endiandata[8], ((uint32_t*)pdata)[8] );
                be32enc( &endiandata[9], ((uint32_t*)pdata)[9] );
                be32enc( &endiandata[10], ((uint32_t*)pdata)[10] );
                be32enc( &endiandata[11], ((uint32_t*)pdata)[11] );
                be32enc( &endiandata[12], ((uint32_t*)pdata)[12] );
                be32enc( &endiandata[13], ((uint32_t*)pdata)[13] );
                be32enc( &endiandata[14], ((uint32_t*)pdata)[14] );
                be32enc( &endiandata[15], ((uint32_t*)pdata)[15] );
                be32enc( &endiandata[16], ((uint32_t*)pdata)[16] );
                be32enc( &endiandata[17], ((uint32_t*)pdata)[17] );
                be32enc( &endiandata[18], ((uint32_t*)pdata)[18] );
                be32enc( &endiandata[19], ((uint32_t*)pdata)[19] );
                be32enc( &endiandata[20], ((uint32_t*)pdata)[20] );
                be32enc( &endiandata[21], ((uint32_t*)pdata)[21] );
                be32enc( &endiandata[22], ((uint32_t*)pdata)[22] );
                be32enc( &endiandata[23], ((uint32_t*)pdata)[23] );
                be32enc( &endiandata[24], ((uint32_t*)pdata)[24] );
                be32enc( &endiandata[25], ((uint32_t*)pdata)[25] );
                be32enc( &endiandata[26], ((uint32_t*)pdata)[26] );
                be32enc( &endiandata[27], ((uint32_t*)pdata)[27] );
                be32enc( &endiandata[28], ((uint32_t*)pdata)[28] );
                be32enc( &endiandata[29], ((uint32_t*)pdata)[29] );
                be32enc( &endiandata[30], ((uint32_t*)pdata)[30] );
                be32enc( &endiandata[31], ((uint32_t*)pdata)[31] );


//	for (int kk=0; kk < 32; kk++) {
//		be32enc(&endiandata[kk], ((uint32_t*)pdata)[kk]);
//	};
#ifdef DEBUG_ALGO
	printf("[%d] Htarg=%X\n", thr_id, Htarg);
#endif
   for (int m=0; m < 6; m++) {
      if (Htarg <= htmax[m]) {
        uint32_t mask = masks[m];
	do {
	   pdata[19] = ++n;
	   be32enc(&endiandata[19], n);
	   x13hash(hash64, endiandata);
#ifndef DEBUG_ALGO
	   if (!(hash64[7] & mask))
           { 
              if ( fulltest(hash64, ptarget) )
              {
	     	*hashes_done = n - first_nonce + 1;
	 	return true;
              }
//                                   else
//                                  {
//                                      applog(LOG_INFO, "Result does not validate on CPU!");
//                                  }
            }
                                   
#else
	    if (!(n % 0x1000) && !thr_id) printf(".");
		if (!(hash64[7] & mask)) {
			printf("[%d]",thr_id);
			if (fulltest(hash64, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				return true;
			}
		}
#endif
	} while (n < max_nonce && !work_restart[thr_id].restart);
			// see blake.c if else to understand the loop on htmax => mask
	break;
     }
  }

  *hashes_done = n - first_nonce + 1;
  pdata[19] = n;
  return 0;
}


bool register_x13_algo( algo_gate_t* gate )
{
  gate->aes_ni_optimized = true;
  init_x13_ctx();
  gate->scanhash = (void*)&scanhash_x13;
  gate->hash     = (void*)&x13hash;
  gate->hash_alt = (void*)&x13hash_alt;
  gate->get_max64 = (void*)&get_max64_0x3ffff;
  return true;
};

