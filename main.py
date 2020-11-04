# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import rlwe_he_scheme_updated as rlwe_updated
import numpy as np

if __name__ == '__main__':
    # Scheme's parameters
    # polynomial modulus degree
    n = 2 ** 2
    # ciphertext modulus
    q = 2 ** 14
    # plaintext modulus
    t = 2
    # base for relin_v1
    T = int(np.sqrt(q)) 
    #modulusswitching modulus
    p = q ** 3

    # polynomial modulus
    poly_mod = np.array([1] + [0] * (n - 1) + [1])
    
    #standard deviation for the error in the encryption
    std1 = 1
    #standard deviation for the error in the evaluateKeyGen_v2
    std2 = 1

    # Keygen
    pk, sk = rlwe_updated.keygen(n, q, poly_mod, std1)

    #EvaluateKeygen_version1
    rlk0_v1, rlk1_v1 = rlwe_updated.evaluate_keygen_v1(sk, n, q, T, poly_mod, std1)

    #EvaluateKeygen_version2
    rlk0_v2, rlk1_v2 = rlwe_updated.evaluate_keygen_v2(sk, n, q, poly_mod, p, std2)
 
    # Encryption
    pt1, pt2 = [1, 0, 1, 1], [1, 1, 0, 1]
    cst1, cst2 = [0, 1, 1, 0], [0, 1, 0, 0]

    ct1 = rlwe_updated.encrypt(pk, n, q, t, poly_mod, pt1, std1)
    ct2 = rlwe_updated.encrypt(pk, n, q, t, poly_mod, pt2, std1)

    print("[+] Ciphertext ct1({}):".format(pt1))
    print("")
    print("\t ct1_0:", ct1[0])
    print("\t ct1_1:", ct1[1])
    print("")
    print("[+] Ciphertext ct2({}):".format(pt2))
    print("")
    print("\t ct1_0:", ct2[0])
    print("\t ct1_1:", ct2[1])
    print("")

    # Evaluation
    ct3 = rlwe_updated.add_plain(ct1, cst1, q, t, poly_mod)
    ct4 = rlwe_updated.mul_plain(ct2, cst2, q, t, poly_mod)
    #ct5 = (ct1 + cst1) + (cst2 * ct2)
    ct5 = rlwe_updated.add_cipher(ct3, ct4, q, poly_mod)
    # ct6 = ct1 * ct2
    ct6 = rlwe_updated.mul_cipher_v1(ct1, ct2, q, t, T, poly_mod, rlk0_v1, rlk1_v1)
    ct7 = rlwe_updated.mul_cipher_v2(ct1, ct2, q, t, p, poly_mod, rlk0_v2, rlk1_v2)
    # Decryption
    decrypted_ct3 = rlwe_updated.decrypt(sk, q, t, poly_mod, ct3)
    decrypted_ct4 = rlwe_updated.decrypt(sk, q, t, poly_mod, ct4)
    decrypted_ct5 = rlwe_updated.decrypt(sk, q, t, poly_mod, ct5)
    decrypted_ct6 = rlwe_updated.decrypt(sk, q, t, poly_mod, ct6)
    decrypted_ct7 = rlwe_updated.decrypt(sk, q, t, poly_mod, ct7)
    
    print("[+] Decrypted ct3=(ct1 + {}): {}".format(cst1, decrypted_ct3))
    print("[+] Decrypted ct4=(ct2 * {}): {}".format(cst2, decrypted_ct4))
    print("[+] Decrypted ct5=(ct1 + {} + {} * ct2): {}".format(cst1, cst2, decrypted_ct5))
    print("[+] pt1 + {} + {} * pt2): {}".format(cst1, cst2, rlwe_updated.polyadd(
                                                rlwe_updated.polyadd(pt1, cst1, t, poly_mod),
                                                rlwe_updated.polymul(cst2, pt2, t ,poly_mod),
                                                t, poly_mod)))
    print("[+] Decrypted ct6=(ct1 * ct2): {}".format(decrypted_ct6))
    print("[+] Decrypted ct7=(ct1 * ct2): {}".format(decrypted_ct7))
    print("[+] pt1 * pt2: {}".format(rlwe_updated.polymul(pt1, pt2, t, poly_mod)))
