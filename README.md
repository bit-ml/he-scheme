We implemented in Python a basic homomorphic encryption scheme inspired from [[FV12]](https://eprint.iacr.org/2012/144.pdf).

## Motivation
The starting point of our implementation is this [github gist](https://gist.github.com/youben11/f00bc95c5dde5e11218f14f7110ad289). The motivation behind our code was for us to understand in detail the two techniques of [[FV12]](https://eprint.iacr.org/2012/144.pdf) used for ciphertext multiplication, namely *relinearization* and *modulus-switching.* This essential operation of ciphertext multiplication was missing in the previous implementation. We thought we might share this understanding through a [blog post](https://bit-ml.github.io/blog/post/homomorphic-encryption-toy-implementation-in-python/) as well since it may be of interest to anyone using the [FV12] scheme in [TenSeal](https://github.com/OpenMined/TenSEAL) or [Seal](https://github.com/Microsoft/SEAL) libraries.

## Disclaimer
Our toy implementation is not meant to be secure or optimized for efficiency. We did it to better understand the inner workings of the [FV12] scheme, so you can use it as a learning tool.

## How to run
The rlwe_he_scheme_updated.py file contains the algorithms of the HE scheme. You can run main.py to play with computing on encrypted data. Have fun! :smile:
