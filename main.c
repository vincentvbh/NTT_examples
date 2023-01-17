
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "NTT_params.h"



#include "gen_table.h"
#include "tools.h"
#include "naive_mult.h"
#include "ntt_c.h"


int16_t table[2048];

int main(void){


    int16_t a[ARRAY_N];
    int16_t b[ARRAY_N];
    int16_t ref[ARRAY_N];
    int16_t res[ARRAY_N];

    int16_t mod, omega, scale, twiddle, t;

    struct compress_profile profile;
    profile.compressed_layers = 3;
    profile.merged_layers[0] = 3;
    profile.merged_layers[1] = 3;
    profile.merged_layers[2] = 2;

    mod = Q1;
    for(size_t i = 0; i < ARRAY_N; i++){
        t = rand();
        cmod_int16(a + i, &t, &mod);
        t = rand();
        cmod_int16(b + i, &t, &mod);
    }

    twiddle = -1;
    mod = Q1;
    naive_mulR(ref,
        a, b,
	ARRAY_N, &twiddle,
	&mod,
	sizeof(int16_t),
	addmod_int16,
	mulmod_int16
    );

    omega = omegaQ1;
    scale = 1;
    mod = Q1;
    gen_twist_table_generic(table,
        &scale, &omega,
	&mod,
	sizeof(int16_t),
	mulmod_int16
    );

    mod = Q1;
    point_mul(a,
	a, table,
	ARRAY_N, 1,
	&mod,
	sizeof(int16_t),
	mulmod_int16
    );
    
    mod = Q1;
    point_mul(b,
	b, table,
	ARRAY_N, 1,
	&mod,
	sizeof(int16_t),
	mulmod_int16
    );


    scale = 1;
    mod = Q1;
    twiddle = omegaQ1;
    mulmod_int16(&omega, &twiddle, &twiddle, &mod);
    gen_streamlined_CT_table_generic(table,
        &scale, &omega,
	&mod,
	sizeof(int16_t),
	mulmod_int16,
	&profile, 0
    );

    mod = Q1;
    // TODO: change the initial 3-layer merge into simple C code
    compressed_CT_NTT_generic(a,
        0, 0,
	table,
	&mod,
	&profile,
	sizeof(int16_t),
	m_layer_CT_butterfly_int16
    );
    compressed_CT_NTT_generic(a,
        1, 2,
	table,
	&mod,
	&profile,
	sizeof(int16_t),
	m_layer_CT_butterfly_int16
    );
    compressed_CT_NTT_generic(b,
        0, 2,
        table,
        &mod,
        &profile,
        sizeof(int16_t),
        m_layer_CT_butterfly_int16
    );

    mod = Q1;
    point_mul(res,
        a, b,
	ARRAY_N, 1,
        &mod,
	sizeof(int16_t),
        mulmod_int16
    );

    scale = 1;
    mod = Q1;
    twiddle = invomegaQ1;
    mulmod_int16(&omega, &twiddle, &twiddle, &mod);
    gen_streamlined_inv_CT_table_generic(table,
        &scale, &omega,
        &mod,
        sizeof(int16_t),
        mulmod_int16,
	expmod_int16,
	&profile, 0
    );

    compressed_CT_iNTT_generic(res,
        0, 2,
        table,
        &mod,
        &profile,
        sizeof(int16_t),
        m_layer_CT_ibutterfly_int16
    );

    scale = invNQ1;
    mod = Q1;
    for(size_t i = 0; i < ARRAY_N; i++){
       mulmod_int16(res + i, res + i, &scale, &mod);
    }

    omega = invomegaQ1;
    scale = 1;
    mod = Q1;
    gen_twist_table_generic(table,
        &scale, &omega,
	&mod,
	sizeof(int16_t),
	mulmod_int16
    );

    mod = Q1;
    point_mul(res,
	res, table,
	ARRAY_N, 1,
	&mod,
	sizeof(int16_t),
	mulmod_int16
    );


    for(size_t i = 0; i < ARRAY_N; i++){
        assert(ref[i] == res[i]);
        //printf("%4zu: %12d, %12d\n", i, ref[i], res[i]);
    } 

    printf("Test finished!\n");


    return 0;

}



