//  racc_api.c
//  Copyright (c) 2023 Raccoon Signature Team. See LICENSE.

//  === Raccoon signature scheme -- NIST KAT Generator API.

#include "api.h"
#include "racc_core.h"
#include "racc_serial.h"
#include "xof_sample.h"
#include <string.h>

//  Generates a keypair - pk is the public key and sk is the secret key.

int	crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
	racc_pk_t r_pk; //  internal-format public key
	racc_sk_t r_sk; //  internal-format secret key
	racc_core_keygen(&r_pk, &r_sk); //  generate keypair
	//  serialize
	if (CRYPTO_PUBLICKEYBYTES != racc_encode_pk(pk, &r_pk)
		|| CRYPTO_SECRETKEYBYTES != racc_encode_sk(sk, &r_sk))
		return (-1);
	return (0);
}

//  Sign a message: sm is the signed message, m is the original message,
//  and sk is the secret key.

int	crypto_sign(unsigned char *sm, unsigned long long *smlen,
		const unsigned char *m, unsigned long long mlen,
		const unsigned char *sk)
{
	uint8_t		mu[RACC_MU_SZ];
	uint64_t	cmul;
	size_t		sig_sz;

	racc_sk_t r_sk;   //  internal-format secret key
	racc_sig_t r_sig; //  internal-format signature
	cmul = 0;
	//  deserialize secret key
	if (CRYPTO_SECRETKEYBYTES != racc_decode_sk(&r_sk, sk))
		return (-1);
	xof_chal_mu(mu, r_sk.pk.tr, m, mlen); //  compute mu
	//  several trials may be needed in case of signature size overflow
	do
	{
    #ifdef BATCH_MODE
		cmul += racc_core_batch_sign(&r_sig, mu, &r_sk); 
	#else
    	cmul += racc_core_sign(&r_sig, mu, &r_sk);
	#endif
		
		//cmul += racc_core_batch_sign(&r_sig, mu, &r_sk); //  create signature
		//  The NIST API expects an "envelope" consisting of the message
		//  together with signature. we put the signature first.
		sig_sz = racc_encode_sig(sm, CRYPTO_BYTES, &r_sig);
	} while (sig_sz == 0);
	memset(sm + sig_sz, 0, CRYPTO_BYTES - sig_sz); //  zero padding
	memcpy(sm + CRYPTO_BYTES, m, mlen);            //  add the message
	*smlen = mlen + CRYPTO_BYTES;
	return (cmul);
}

//  Verify a message signature: m is the original message, sm is the signed
//  message, pk is the public key.

int	crypto_sign_open(unsigned char *m, unsigned long long *mlen,
		const unsigned char *sm, unsigned long long smlen,
		const unsigned char *pk)
{
	size_t	m_sz;
	uint8_t	mu[RACC_MU_SZ];

	racc_pk_t r_pk;   //  internal-format public key
	racc_sig_t r_sig; //  internal-format signature
	//  deserialize public key, signature with a consistency check
	if (smlen < CRYPTO_BYTES || CRYPTO_PUBLICKEYBYTES != racc_decode_pk(&r_pk,
			pk) || CRYPTO_BYTES != racc_decode_sig(&r_sig, sm))
		return (-1);
	m_sz = smlen - CRYPTO_BYTES;
	//  compute mu
	xof_chal_mu(mu, r_pk.tr, sm + CRYPTO_BYTES, m_sz);
	//  verification
	if (!racc_core_verify(&r_sig, mu, &r_pk))
		return (-1);
	//  store the length and move the "opened" message
	memcpy(m, sm + CRYPTO_BYTES, m_sz);
	*mlen = m_sz;
	return (0);
}
