# pyNTRUEncrypt
Python NTRUEncrypt implementation.

It is licensed under GPL v2 because that is what the existing patents are open-sourced under, assuming my bare implementation infringes.

Until I get the chance to refactor it, it's all in one file and exposes much more than it needs to.

Of note, `NTRUEncrypt` takes an `NTRUParams` object, a `ConvModQ` object (which is a public key parameter), and a message `m`, which is assumed to be a string.

The real magic takes place in `NTRUBlockEncrypt` and `NTRUBlockDecrypt` which take the same parameters, but assumes `m` is already encoded. Right now, we convert the string into binary and split it among however many convolution polynomials are required. In the future, I would like to encode this as a sequence of trinary polynomials because `m` is  `mod p`, where `p = 3` in the standard implemenation.

Currently, an NTRU key can be automatically generated for you (defaulting to 256-bit) or can be explicitly declared by passing an `NTRUParams` and two polynomials, `f` and `g`. If passed, they will be tried but will be replaced if found to not possess an inverse properly.

Rings can be generated using `NTRUParams`. Current choices are 112-bit, 128-bit, 192-bit, and 256-bit. Additionally, you may choose to emphasize space, speed, or a combination. These parameters are from the latest (to my knowledge, at time of development) whitepaper produced by the NTRU folks, found here: https://www.securityinnovation.com/uploads/Crypto/params.pdf.

# Should I use this?

Oh, God no! This was purely an exercise in mathematical logic as applied to cryptosystems. While easily ported to other languages, this is a naive implemenation and is bound to be open to numerous side-channel attacks.

# How fast is it?

You know, I was going to say it was really slow, but it turns out that it's only slow under CPython. To encrypt and then decrypt `"Hello World!"` is CPython, it takes about `0.84s`. In PyPy, it drops down to `0.077s` which is actually somewhat useable.


