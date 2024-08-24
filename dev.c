#include <stdio.h>
#include <pbc/pbc.h>
#include <gmp.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define num_user 1    //Define the global variable for total number of users

typedef struct Node {
    mpz_t name; //Index number or Identity of node
    mpf_t secret_value; // Arbitrary-precision float value according to lagrange's algorithm
    long unsigned int threshold; // Presenting the total number of initial dummy nodes
} Node;

typedef struct public {
    pbc_param_t param; // pairing parameter
    element_t G1; // g1 for group G1
    element_t G2; // g2 for group G2
    mpz_t p; // order of bilinear pairing
    mpz_t k; // security parameter
    mpz_t a; // total number of attributes
    pairing_t pairing; // bilinear pairing according to parameter param
    element_t identity_g1; // identity for group G1
    element_t identity_g2; // identity for group G2
    element_t *generators_G1; // generators of group G1
    element_t *generators_G2; // generators of group G2

    char *message; // public message
    element_t wf; // water functions
    element_t *gatt; // gatt[j] = g1^secret[j] and component of group public key
    element_t omega; // random component of group public key
    element_t *u; // groth sahai proof vector and component of group public key
    element_t *u_dash; // groth sahai proof vector and component of group public key
    element_t *v; // groth sahai proof vector and component of group public key

    element_t vT; // vT = G1^sT
    mpz_t d; // Total number of dummy nodes
    Node *T_d; // Array of dummy nodes depicting the Extension tree
    Node *T_trial; // Temporary tree used for the rough purposes and not part of research paper
    mpz_t max; // Total number of the nodes in the Extension tree
    mpz_t temp1; // Total number of non-dummy nodes in the Extension tree
    mpz_t dummy; // Total number of dummy nodes in the Extension tree
}
public;

public ret_setup;
// The ret_setup is variable for the structure public which means the public values

typedef struct secret {
    mpz_t *S; // secret values of all attributes
    element_t gamma; // random component of issuing key
    element_t alpha1; // Open user key
    element_t alpha1_dash; // attribute tracing key
    Node *T; // Secret valued tree made using the secret values of attributes and dummy nodes
} secret;

secret private;
// The private is variable for the structure secret which means the secret values

typedef struct user_output {
    long unsigned int i; // user's name
    mpz_t upki; // user public key
    mpz_t uski; // user secret key
    element_t Ai; // A security function part of digital signature and membership certificate
    element_t Xi_1; // Xi_1 = G1 ^ xi and part of membership certificate
    element_t Xi_2; // Xi_1 = G2 ^ xi and part of membership certificate
    element_t yi_dash_dash; // A random value chose by user

    mpz_t signature; // digital signature used for verification of user
    mpz_t rsa1; // RSA encryption key belong to user
    mpz_t rsa2; // Public RSA encryption key
    element_t yi; // A random value chose by the user and part of membership certificate
    element_t Yi; // Yi = G1 ^ yi
    element_t *Ti; // Attribute certificate of the user
    element_t xi; // A random value chose by the user
    long unsigned int a; // Total number of attributes assigned to the user
    long unsigned int *A; // Array of attributes assigned to user

    Node *T; // Secret valued tree Array made using the secret values of user specific attributes
    Node *T_d; // Secret valued tree Array made using the secret values of user specific dummy nodes

    element_t vT; // vT = G1 * sT specifically derived by the user

    element_t rho1; // equal to Ai above
    element_t rho2; // equal to yi above
    element_t *rho3; // equal to Xi_1 and Xi_2
    element_t rho4; // part of signature
    element_t rho5; // rho5 = h ^ sT2 (i.e. h is the first generator of group G1)
    element_t rho6; // ID of user
    element_t *rho7; // rho7[0] = G1 ^ z and rho7[1] = G2 ^ z
    element_t rho8; // part of signature
    element_t rho9; // rho9 = G2 ^ r
    element_t *sigma1; // commitment of rho1
    element_t *sigma2_1; // commitment of rho2[0]
    element_t *sigma2_2; // commitment of rho2[1]
    element_t *sigma3_1; // commitment of rho3[0]
    element_t *sigma3_2; // commitment of rho3[1]
    element_t *sigma4; // commitment of rho4
    element_t *sigma5; // commitment of rho5
} user_output;

// The user_output structure is used to store data which is at user end

user_output U[num_user];
// The U is the structure user_output

typedef struct reg {
    long unsigned int i; // user's register name
    mpz_t upki; // user public key
    long unsigned int a; // Total number of attributes assigned to the user
    long unsigned int *A; // Array of attributes assigned to user
    element_t Ai; // A security function part of digital signature and membership certificate
    element_t Xi_1; // Xi_1 = G1 ^ xi and part of membership certificate
    element_t Xi_2; // Xi_1 = G2 ^ xi and part of membership certificate
    element_t Yi; // Yi = G1 ^ yi
    mpz_t signature; // digital signature used for verification of user
    mpz_t rsa2; // Public RSA encryption key
    mpz_t rsa3; // RSA encryption key belong to group manager

    mpz_t sT; // secret value of the Extension Tree
    mpz_t sT1; // Sum of all secret values of attributes * lagrange's constant
    mpz_t sT2; // Sum of all secret values of dummy nodes * lagrange's constant
    element_t r; // Random element for encryption
    element_t z; // Random element for encryption
} reg;

// The reg structure is used to store data which is at Register/user

reg R[num_user];
// The R is the structure reg as (index number + 1) as user's identity

// random_prime_bits() function to make a prime number using given number of bites
void random_prime_bits(mpz_t result, mpz_t n) {
    // Declare variable random
    mpz_t random; //REF-01
    mpz_init(random); //REF-02
    // Declare state
    gmp_randstate_t state; //REF-04
    // Initialize state
    gmp_randinit_default(state); //REF-05
    // Set initial seed value to state. We pass seed as a random number defined as (random_number + 1)*(another-random-number)
    // + 1. 1 was added to avoid seed = 0.
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1)); //REF-06
    // If number of bits <= 1, then no prime exists for that number of bits
    if (mpz_cmp_ui(n, 1) <= 0) {
        //REF-07
        printf("NO PRIME EXISTS\n");
    } else {
        mpz_t lower_limit; //REF-01
        mpz_init(lower_limit); //REF-02
        mpz_ui_pow_ui(lower_limit, 2, mpz_get_ui(n) - 1); //REF-08  REF-09
        // Loop till we find a prime number of n bits
        while (1) {
            // Store a random value from 0 to (2^n)-1 in the variable random
            mpz_urandomb(random, state, mpz_get_ui(n)); //REF-10  REF-09
            // Then check whether random number is prime we use the probabilistic function
            if (mpz_cmp(random, lower_limit) > 0 && mpz_probab_prime_p(random, mpz_get_ui(n))) {
                //REF-11  REF-12  //REF-09
                // If random is of n-bits and is prime (probably or for sure), then set result to the random number generated and return
                mpz_set(result, random); //REF-13
                break;
            }
        }
    }
}

// generate_rsa_keys() function is used to generate keys for applying the RSA crypto algorithm
void generate_rsa_keys(mpz_t n, mpz_t d, mpz_t e, gmp_randstate_t state) {
    mpz_t p, q, phi; //REF-01
    mpz_inits(p, q, phi, NULL); //REF-02
    mpz_set(e, ret_setup.p); // Common public exponent       REF-13

    // Generate two large primes p and q
    mpz_urandomb(p, state, mpz_get_ui(ret_setup.k)); //REF-10
    mpz_urandomb(q, state, mpz_get_ui(ret_setup.k)); //REF-10
    mpz_nextprime(p, p); //REF-41
    mpz_nextprime(q, q); //REF-41

    // Calculate n = p * q
    mpz_mul(n, p, q); //REF-42
    // Calculate phi = (p-1)*(q-1)
    mpz_sub_ui(p, p, 1); //REF-43
    mpz_sub_ui(q, q, 1); //REF-43
    mpz_mul(phi, p, q); //REF-42
    // Calculate d = e^(-1) mod phi
    mpz_invert(d, e, phi); //REF-44
    // Clearing the temporary memory
    mpz_clears(p, q, phi, NULL);
}

// insertionSort() function is the sorting algorithm using insertion sort for proper orientation for storing attributes
void insertionSort(long unsigned int arr[], long unsigned int n) {
    int key, j;
    for (int i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;

        // Move elements of arr[0..i-1] that are greater than key to one position ahead of their current position
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

// initializeNode() function is used to intialize the node for avoiding segmentation faults and using node structure easily
void initializeNode(Node *node) {
    mpz_init(node->name); //REF-02
    mpf_init(node->secret_value); //REF-50
    node->threshold = 1; // Set default threshold = 1
}

void tree_making(mpz_t dummy, mpz_t temp1, mpz_t max, long unsigned int A[], Node *T_d, mpz_t *S, long unsigned int z) {
    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node)); //REF-09
    long unsigned int i = 0;
    for (i = 0; i < mpz_get_ui(temp1); i++) {
        //REF-09
        initializeNode(&T[i]);
    }
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        //REF-09
        i = j - mpz_get_ui(dummy); //REF-09
        mpf_set_z(T[j].secret_value, S[i]); //REF-48
    }
    mpz_set_ui(T[0].name, 1); //REF-03
    long unsigned int k = 1, m = 2, l = 0;
    while (m <= mpz_get_ui(max)) {
        //REF-09
        mpz_set_ui(T[k].name, m); //REF-03
        mpz_set_ui(T[k + 1].name, m + 1); //REF-03
        mpz_set_ui(T_d[l].name, m + 2); //REF-03
        k = k + 2;
        l++;
        m = m + 3;
    }
    i = 0;
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        //REF-09
        if (j - mpz_get_ui(dummy) == A[i]) {
            //REF-09
            i++;
        } else {
            initializeNode(&T[j]);
        }
    }
    i = mpz_get_ui(temp1) - 1; //REF-09
    printf("\n");
    mpf_t temp2, temp3, temp4; //REF-49
    mpf_inits(temp2, temp3, temp4, NULL); //REF-50
    mpz_t temp5; //REF-01
    mpz_init(temp5); //REF-02
    while (i > 0) {
        if (mpz_cmp_ui(T[i].name, 0) == 0 && mpz_cmp_ui(T[i - 1].name, 0) == 0) {
            //REF-07
            initializeNode(&T[i / 2 - 1]);
            initializeNode(&T_d[i / 2 - 1]);
        } else if (mpz_cmp_ui(T[i - 1].name, 0) == 0) {
            //REF-07
            mpf_set_z(temp2, T_d[i / 2 - 1].name); //REF-48
            mpf_mul(temp2, temp2, T[i].secret_value); //REF-51
            mpf_set_z(temp3, T[i].name); //REF-48
            mpf_mul(temp3, temp3, T_d[i / 2 - 1].secret_value); //REF-51
            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3); //REF-52
        } else if (mpz_cmp_ui(T[i].name, 0) == 0) {
            //REF-03
            mpf_set_z(temp2, T_d[i / 2 - 1].name); //REF-48
            mpf_div_ui(temp2, temp2, 2); //REF-59
            mpf_set_z(temp3, T[i - 1].name); //REF-48
            mpf_div_ui(temp3, temp3, 2); //REF-59
            mpf_mul(temp2, T[i - 1].secret_value, temp2); //REF-51
            mpf_mul(temp3, T_d[i / 2 - 1].secret_value, temp3); //REF-51
            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3); //REF-52
        } else {
            mpz_mul(temp5, T[i - 1].name, T[i].name); //REF-42
            mpf_set_z(temp2, temp5); //REF-48
            mpf_mul(temp2, T_d[i / 2 - 1].secret_value, temp2); //REF-51

            mpz_mul(temp5, T[i].name, T_d[i / 2 - 1].name); //REF-42
            mpf_set_z(temp3, temp5); //REF-48
            mpf_mul(temp4, temp3, T[i - 1].secret_value); //REF-51

            mpf_add(temp2, temp2, temp4); //REF-60
            mpf_div_ui(temp2, temp2, 2); //REF-59

            mpz_mul(temp5, T_d[i / 2 - 1].name, T[i - 1].name); //REF-42
            mpf_set_z(temp3, temp5); //REF-48
            mpf_mul(temp3, temp3, T[i].secret_value); //REF-51

            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3); //REF-52
        }
        i = i - 2;
    }
    if (z == 0)
        ret_setup.T_trial = T;
    else
        U[z - 1].T = T;

    i = mpz_get_ui(temp1) - 1; //REF-09
    while (i > 0) {
        long unsigned int j = (i / 2) - 1;
        mpf_mul_ui(temp4, private.T[i].secret_value, 2); //REF-53
        mpf_sub(T_d[j].secret_value, temp4, private.T[i - 1].secret_value); //REF-52
        i = i - 2;
    }
    m = 2, l = 0;
    while (m <= mpz_get_ui(max)) {
        //REF-09
        mpz_set_ui(T_d[l].name, m + 2); //REF-03
        l++;
        m = m + 3;
    }
}

void verify_tree(long unsigned int A[], mpz_t dummy, mpz_t temp1, mpz_t max, Node *T_d, unsigned long int z) {
    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node));
    mpz_t *S = (mpz_t *) malloc(mpz_get_ui(ret_setup.a) * sizeof(mpz_t)); //REF-01
    long unsigned int i = 0;
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        //REF-09
        long unsigned int k = j - mpz_get_ui(dummy); //REF-09
        if (k == A[i]) {
            i++;
            mpz_init_set(S[k], private.S[k]); //REF-54
        } else
            mpz_init(S[k]); //REF-02
    }
    tree_making(dummy, temp1, max, A, T_d, S, z);
    if (z == 0)
        T = ret_setup.T_trial;
    else
        T = U[z - 1].T;
    gmp_printf("Root  : %.Ff\n", T[0].secret_value); //REF-14
    long unsigned int l = 0, m = 1, k = 1;
    while (m < mpz_get_ui(max)) {
        //REF-09
        gmp_printf("Node %Zd: %.Ff\n", T[k].name, T[k].secret_value); //REF-14
        gmp_printf("Node %Zd: %.Ff\n", T[k + 1].name, T[k + 1].secret_value); //REF-14
        k = k + 2;
        gmp_printf("Node %Zd: %.Ff\n", T_d[l].name, T_d[l].secret_value); //REF-14
        l++;
        m = m + 3;
    }
    mpz_t temp5; //REF-01
    mpz_init(temp5); //REF-02
    element_t groot, vT, tem1, tem2, sT; //REF-19
    element_init_G1(groot, ret_setup.pairing); //REF-20
    element_init_G2(vT, ret_setup.pairing); //REF-21
    element_init_GT(tem1, ret_setup.pairing); //REF-22
    element_init_GT(tem2, ret_setup.pairing); //REF-22
    element_init_Zr(sT, ret_setup.pairing); //REF-31
    mpz_set_f(temp5, T[0].secret_value); //REF-55
    mpz_abs(temp5, temp5); //REF-56
    element_set_mpz(sT, temp5); //REF-57
    element_pow_zn(groot, ret_setup.G1, sT); //REF-34
    element_set(vT, ret_setup.vT); //REF-28
    element_pairing(tem1, groot, ret_setup.G2); //REF-38
    element_pairing(tem2, ret_setup.G1, vT); //REF-38
    if (element_cmp(tem1, tem2) == 0) {
        //REF-27
        if (z != 0) {
            printf("\n-----Signing verification valid-----\n");
        } else
            printf("\n-----Buildtree algorithm valid-----\n");
    } else
        printf("\nx-x-x-Buildtree algorithm invalid-x-x-x\n");

    Node *T_dtrial = (Node *) malloc(mpz_get_ui(dummy) * sizeof(Node));
    for (i = 0; i < mpz_get_ui(dummy); i++) {
        //REF-09
        initializeNode(&T_dtrial[i]);
    }
    tree_making(dummy, temp1, max, A, T_dtrial, S, 0);

    mpz_t sT1, sT2, s; //REF-01
    mpz_inits(sT1, sT2, s, NULL); //REF-02
    mpz_set_f(sT1, ret_setup.T_trial[0].secret_value); //REF-55
    for (i = 0; i < mpz_get_ui(ret_setup.a); i++) {
        //REF-09
        mpz_set_ui(S[i], 0); //REF-03
    }

    tree_making(dummy, temp1, max, A, T_d, S, 0);
    mpz_set_f(sT2, ret_setup.T_trial[0].secret_value); //REF-55
    mpz_add(s, sT1, sT2); //REF-40
    if (z != 0) {
        mpz_inits(R[z - 1].sT, R[z - 1].sT1, R[z - 1].sT2, NULL); //REF-02
        if (mpz_cmp_ui(s, 0) < 0) {
            //REF-07
            mpz_abs(R[z - 1].sT, s); //REF-56
            if (mpz_cmp_ui(sT1, 0) < 0) {
                //REF-07
                mpz_abs(R[z - 1].sT1, sT1); //REF-56
                mpz_neg(R[z - 1].sT2, sT2); //REF-61
            } else {
                mpz_abs(R[z - 1].sT2, sT2); //REF-56
                mpz_neg(R[z - 1].sT1, sT1); //REF-61
            }
        } else {
            mpz_set(R[z - 1].sT1, sT1); //REF-13
            mpz_set(R[z - 1].sT, s); //REF-13
            mpz_set(R[z - 1].sT2, sT2); //REF-13
        }
    }
    gmp_printf("s  : %Zd\n", s); //REF-14
    gmp_printf("st1  : %Zd\n", sT1); //REF-14
    gmp_printf("st2  : %Zd\n", sT2); //REF-14
}

// setup() function is used to generate parameters and make the setup algorithm
void setup(public *retval, mpz_t k) {
    printf("\n----------------Setup Algorithm------------------\n");
    mpz_t p; //REF-01
    mpz_init(p); //REF-02
    mpz_init(retval->a); //REF-02
    // Total Attributes = security parameter
    mpz_set(retval->a, k); //REF-13
    // generating prime order of bilinear pairing
    random_prime_bits(p, k);
    // Printing the bilinear group prime order and the total number of attribute
    gmp_printf("P(order of bilinear group) = %Zd\n", p); //REF-14
    gmp_printf("Total Attributes = %Zd\n", retval->a); //REF-14

    // Intializing pairing for security parameters
    pairing_t pairing; //REF-15
    pbc_param_t param; //REF-16
    pbc_param_init_a1_gen(param, p); //REF-17
    pairing_init_pbc_param(pairing, param); //REF-18

    // setting various entities as per reasearch paper
    element_t g1, g2, gt1, identity_g1, identity_g2, temp_g1, temp_g2; //REF-19
    element_init_G1(g1, pairing); //REF-20
    element_init_G2(g2, pairing); //REF-21
    element_init_GT(gt1, pairing); //REF-22
    element_init_G1(temp_g1, pairing); //REF-20
    element_init_G2(temp_g2, pairing); //REF-21
    element_init_G1(identity_g1, pairing); //REF-20
    element_init_G2(identity_g2, pairing); //REF-21

    element_set0(identity_g2); //REF-23
    element_set0(identity_g1); //REF-23

    mpz_t required, gen; //REF-01
    mpz_inits(required, gen, NULL); //REF-02
    mpz_add_ui(required, k, 2); //REF-24

    // storing the generators of the group G1 and G2
    element_t *generators_g1 = (element_t *) malloc(sizeof(element_t) * (mpz_get_ui(required))); //REF-19
    element_t *generators_g2 = (element_t *) malloc(sizeof(element_t) * (mpz_get_ui(required))); //REF-19

    for (unsigned long int i = 0; i < mpz_get_ui(required); i++) {
        element_init_G1(generators_g1[i], pairing); //REF-20
        element_init_G2(generators_g2[i], pairing); //REF-21
    }
    unsigned long long int index = 0;
    do {
        element_random(g1); //REF-25
        element_pow_mpz(temp_g1, g1, p); //REF-26
        if (element_cmp(temp_g1, identity_g1) == 0) {
            //REF-27
            mpz_add_ui(gen, gen, 1); //REF-23
            element_set(generators_g1[index], g1); //REF-28
            index++;
        }
    } while (mpz_cmp(gen, required)); //REF-11
    index = 0;

    mpz_set_ui(gen, 0); //REF-09
    do {
        element_random(g2); //REF-25
        element_pow_mpz(temp_g2, g2, p); //REF-26
        if (element_cmp(temp_g2, identity_g2) == 0) {
            //REF-27
            mpz_add_ui(gen, gen, 1); //REF-23
            element_set(generators_g2[index], g2); //REF-28
            index++;
        }
    } while (mpz_cmp(gen, required)); //REF-11

    // Storing the necessary things to ret_setup structure for further use
    pbc_param_init_a1_gen(retval->param, p); //REF-17
    pairing_init_pbc_param(retval->pairing, param); //REF-18
    mpz_init(retval->k); //REF-02
    mpz_set(retval->k, k); //REF-13
    mpz_init(retval->p); //REF-02
    mpz_set_ui(retval->p, mpz_get_ui(p)); //REF-03
    element_init_G1(retval->G1, pairing); //REF-20
    element_init_G2(retval->G2, pairing); //REF-21
    element_init_G1(retval->identity_g1, pairing); //REF-20
    element_init_G2(retval->identity_g2, pairing); //REF-21
    element_set0(retval->identity_g1); //REF-23
    element_set0(retval->identity_g2); //REF-23
    retval->generators_G1 = generators_g1;
    retval->generators_G2 = generators_g2;
}

// keygen() function takes an input system parameters params and outputs a group public key gpk, an issuing key ,a user
// opening key ok_user (alpha1) and an attribute tracing key (alpha1_dash)
void keygen(char *M, public *ret_val, secret *sec_val) {
    printf("\n----------------KenGen Algorithm------------------\n");

    gmp_randstate_t state; //REF-04
    gmp_randinit_default(state); //REF-05
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1)); //REF-06

    // Generating the water function
    element_t wf, temp_pw; //REF-19
    element_init_G1(wf, ret_val->pairing); //REF-20
    element_set(wf, ret_val->generators_G1[1]); //REF-28
    element_init_G1(temp_pw, ret_val->pairing); //REF-20
    mpz_t msg_val; //REF-01
    mpz_init(msg_val); //REF-02

    for (unsigned long int i = 0; i < mpz_get_ui(ret_val->k); i++) {
        if (M[i] == '1') {
            mpz_set_ui(msg_val, 1); //REF-03
        } else {
            mpz_set_ui(msg_val, 0); //REF-03
        }
        element_pow_mpz(temp_pw, ret_setup.generators_G1[i + 2], msg_val); //REF-26
        element_mul(wf, wf, temp_pw); //REF-29
    }

    printf("Group Public Key : \n");

    element_printf("  h = %B\n", ret_val->generators_G1[0]); //REF-30

    element_t g1, g2; //REF-19
    element_init_G1(g1, ret_val->pairing); //REF-20
    element_init_G2(g2, ret_val->pairing); //REF-21
    element_random(g1); //REF-25
    element_random(g2); //REF-25
    element_t gamma, omega; //REF-19
    element_init_Zr(gamma, ret_val->pairing); //REF-31
    element_init_G2(omega, ret_val->pairing); //REF-21
    element_random(gamma); //REF-25
    element_mul_zn(omega, g2, gamma); //REF-32

    element_set(ret_val->G1, g1); //REF-28
    element_set(ret_val->G2, g2); //REF-28
    element_init_G2(ret_val->omega, ret_val->pairing); //REF-21
    element_set(ret_val->omega, omega); //REF-28
    element_init_Zr(sec_val->gamma, ret_val->pairing); //REF-31
    element_set(sec_val->gamma, gamma); //REF-28
    element_init_G1(ret_val->wf, ret_val->pairing); //REF-20
    element_set(ret_val->wf, wf); //REF-28
    element_printf("  Omega = %B\nWater Function : %B\n", ret_val->omega, wf); //REF-30
    long unsigned int a = mpz_get_ui(ret_val->a);
    mpz_t *S = (mpz_t *) malloc(sizeof(mpz_t) * a); //REF-01
    element_t *gatt = (element_t *) malloc(sizeof(element_t) * a); //REF-19
    printf("  Generator Attributes : For particlar attributes you can print\n");

    for (unsigned long int i = 0; i < a; i++) {
        mpz_init(S[i]); //REF-02
        while (mpz_cmp_ui(S[i], 0) == 0) //REF-07
            mpz_urandomm(S[i], state, ret_setup.p); //REF-33
        element_init_G1(gatt[i], ret_setup.pairing); //REF-20
        element_pow_mpz(gatt[i], g1, S[i]); //REF-26
    }
    ret_val->gatt = gatt;
    sec_val->S = S;

    element_t alpha1, alpha2, t1, t2, alpha1_dash; //REF-19
    element_init_Zr(alpha1, ret_val->pairing); //REF-31
    element_init_Zr(alpha2, ret_val->pairing); //REF-31
    element_init_Zr(alpha1_dash, ret_val->pairing); //REF-31
    element_init_Zr(t1, ret_val->pairing); //REF-31
    element_init_Zr(t2, ret_val->pairing); //REF-31

    element_random(alpha2); //REF-25
    element_random(alpha1); //REF-25
    element_random(alpha1_dash); //REF-25
    element_random(t1); //REF-25
    element_random(t2); //REF-25

    element_t u1_dash[2], v1[2], u1[2]; //REF-19
    element_t *u = (element_t *) malloc(sizeof(element_t) * 4); //REF-19
    element_t *u_dash = (element_t *) malloc(sizeof(element_t) * 4); //REF-19
    element_t *v = (element_t *) malloc(sizeof(element_t) * 4); //REF-19

    for (int i = 0; i < 2; i++) {
        element_init_G1(u1[i], ret_setup.pairing); //REF-20
        element_init_G1(u1_dash[i], ret_setup.pairing); //REF-20
        element_init_G2(v1[i], ret_setup.pairing); //REF-21
    }

    for (int i = 0; i < 4; i++) {
        element_init_G1(u[i], ret_setup.pairing); //REF-20
        element_init_G1(u_dash[i], ret_setup.pairing); //REF-20
        element_init_G2(v[i], ret_setup.pairing); //REF-21
    }

    element_set(u1[0], g1); //REF-28
    element_pow_zn(u1[1], g1, alpha1); //REF-34
    element_set(u1_dash[0], g1); //REF-28
    element_pow_zn(u1_dash[1], g1, alpha1_dash); //REF-34
    element_set(v1[0], g2); //REF-28
    element_pow_zn(v1[1], g2, alpha2); //REF-34
    element_set(u[0], u1[0]); //REF-28
    element_set(u[1], u1[1]); //REF-28
    element_pow_zn(u[2], u1[0], t1); //REF-34
    element_pow_zn(u[3], u1[1], t1); //REF-34
    element_set(u_dash[0], u1_dash[0]); //REF-28
    element_set(u_dash[1], u1_dash[1]); //REF-28
    element_pow_zn(u_dash[2], u1_dash[0], t1); //REF-34
    element_pow_zn(u_dash[3], u1_dash[1], t1); //REF-34
    element_set(v[0], v1[0]); //REF-28
    element_set(v[1], v1[1]); //REF-28
    element_pow_zn(v[2], v1[0], t2); //REF-34
    element_pow_zn(v[3], v1[1], t2); //REF-34

    ret_setup.u = u;
    ret_setup.v = v;
    ret_setup.u_dash = u_dash;

    printf("  Groth Sahai Proof Vectors :\n");
    for (int i = 0; i < 4; i++) {
        element_printf("\tu[%d]: %B\n", i, u[i]); //REF-30
        element_printf("\tu_dash[%d]: %B\n", i, u_dash[i]); //REF-30
        element_printf("\tv[%d]: %B\n", i, v[i]); //REF-30
    }

    element_printf("Issuing key : \n  Gamma = %B\n"
                   "  Set S = For particular attributes of secrets you can print\n",
                   gamma); //REF-30
    element_printf("Opener key = %B\n", alpha1); //REF-30
    element_printf("Tracer key = %B\n", alpha1_dash); //REF-30
}

void join(public *retval, unsigned long int j, user_output *U_val, reg *rval, secret *sec_val) {
    printf("\n----------------Joint Algorithm for User U%lu------------------\n", j);

    U_val->i = j;
    rval->i = j;

    gmp_randstate_t state; //REF-04
    gmp_randinit_default(state); //REF-05
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1)); //REF-05

    mpz_inits(U_val->upki, rval->upki, U_val->uski, NULL);
    mpz_urandomm(U_val->upki, state, retval->p); //REF-33
    mpz_set(rval->upki, U_val->upki); //REF-13
    mpz_urandomm(U_val->uski, state, retval->p); //REF-28

    element_t yi_dash, g1, Yi_dash; //REF-19
    element_init_Zr(yi_dash, retval->pairing); //REF-31
    element_init_G1(g1, retval->pairing); //REF-20
    element_init_G1(Yi_dash, retval->pairing); //REF-20
    element_random(yi_dash); //REF-25
    element_set(g1, ret_setup.G1); //REF-28
    element_mul_zn(Yi_dash, g1, yi_dash); //REF-32
    element_printf("User U%lu sends :\n  Yi_dash : %B\n", j, Yi_dash); //REF-30

    element_t xi, yi_dash_dash, Yi_dash_dash, Ai, gamma, add; //REF-19
    element_init_Zr(xi, retval->pairing); //REF-31
    element_init_Zr(yi_dash_dash, retval->pairing); //REF-31
    element_init_Zr(U_val->yi_dash_dash, retval->pairing); //REF-31
    element_init_G1(Yi_dash_dash, retval->pairing); //REF-20
    element_init_G1(Ai, retval->pairing); //REF-20
    element_init_Zr(gamma, retval->pairing); //REF-31
    element_init_Zr(add, retval->pairing); //REF-31
    element_random(xi); //REF-25
    element_random(yi_dash_dash); //REF-25
    element_set(U_val->yi_dash_dash, yi_dash_dash); //REF-28
    element_pow_zn(Yi_dash_dash, g1, yi_dash_dash); //REF-34
    element_set(gamma, sec_val->gamma); //REF-28

    element_set(Ai, retval->generators_G1[0]); //REF-28
    element_mul(Ai, Ai, Yi_dash); //REF-29
    element_mul(Ai, Ai, Yi_dash_dash); //REF-29
    element_add(add, gamma, xi); //REF-35
    element_invert(add, add); //REF-36
    element_pow_zn(Ai, Ai, add); //REF-34

    element_t Xi_2; //REF-19
    element_init_G2(Xi_2, retval->pairing); //REF-21
    element_set(Xi_2, retval->G2); //REF-28
    element_pow_zn(Xi_2, Xi_2, xi); //REF-34

    mpz_t a; //REF-01
    mpz_init(a); //REF-02
    while (mpz_get_ui(a) == 0) {
        //REF-09
        mpz_urandomm(a, state, retval->a); //REF-33
    }

    long unsigned int *A = (long unsigned int *) malloc(sizeof(long unsigned int) * mpz_get_ui(a));
    long unsigned int k = 0;
    do {
        int z = 0;
        unsigned long int y = mpz_get_ui(retval->a); //REF-09
        unsigned long int temp1 = rand() % y;
        for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
            //REF-09
            if (temp1 == A[i]) {
                z = 1;
                break;
            }
        }
        if (z == 0) {
            A[k] = temp1;
            k++;
        }
    } while (k < mpz_get_ui(a)); //REF-09

    element_t temp_pow; //REF-19
    element_init_Zr(temp_pow, retval->pairing); //REF-31
    element_t *Ti = (element_t *) malloc(sizeof(element_t) * mpz_get_ui(retval->k)); //REF-19
    for (long unsigned int i = 0; i < mpz_get_ui(retval->a); i++) {
        //REF-09
        element_init_G1(Ti[i], retval->pairing); //REF-20
        element_mul_mpz(temp_pow, add, sec_val->S[i]); //REF-38
        element_pow_zn(Ti[i], retval->generators_G1[0], temp_pow); //REF-34
        element_set0(temp_pow); //REF-23
    }

    element_printf("Group Manager sends :\n"); //REF-30
    element_printf("  y%lu_dash_dash : %B\n", j, yi_dash_dash); //REF-30
    element_printf("  A%lu : %B\n", j, Ai); //REF-30
    element_printf("  X%lu_2 : %B\n", j, Xi_2); //REF-30
    for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
        //REF-09
        long unsigned int x = A[i];
        element_printf("  T%lu_%lu : %B\n", j, x, Ti[x]); //REF-30
    }

    element_t yi, Yi, temp1, temp2, temp3, temp4, temp_g2; //REF-19
    element_init_Zr(yi, retval->pairing); //REF-31
    element_init_G1(Yi, retval->pairing); //REF-20
    element_init_GT(temp1, retval->pairing); //REF-22
    element_init_GT(temp2, retval->pairing); //REF-22
    element_init_GT(temp3, retval->pairing); //REF-22
    element_init_GT(temp4, retval->pairing); //REF-22
    element_init_G2(temp_g2, retval->pairing); //REF-21
    element_add(yi, yi_dash, yi_dash_dash); //REF-35
    element_mul(temp_g2, Xi_2, ret_setup.omega); //REF-29
    element_pairing(temp1, Ai, temp_g2); //REF-38

    element_pairing(temp2, retval->generators_G1[0], retval->G2); //REF-38
    element_pairing(temp3, retval->G1, retval->G2); //REF-38
    element_pow_zn(temp3, temp3, yi); //REF-34
    element_mul(temp4, temp3, temp2); //REF-29

    printf("User U%lu processing : \n", j);
    if (element_cmp(temp1, temp4) == 0) //REF-27
        printf("\n   -----First Verification Done by user-----\n");
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");
    element_pow_zn(Yi, g1, yi); //REF-34

    printf("\n  Digital Signature Sigma_%lu for : \n", j);
    element_printf("\tA%lu : %B\n", j, Ai); //REF-30
    element_printf("\tX%lu_2 : %B\n", j, Xi_2); //REF-30
    element_printf("\tY%lu : %B\n   Processing i.e. DSig_usk%lu :\n", j, Yi, j); //REF-30

    mpz_t hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki; //REF-01
    mpz_inits(hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki, NULL); //REF-02
    element_to_mpz(hash_Ai, Ai); //REF-39
    element_to_mpz(hash_Xi_2, Xi_2); //REF-39
    element_to_mpz(hash_Yi, Yi); //REF-39

    mpz_add(combined_hash, hash_Ai, hash_Xi_2); //REF-40
    mpz_add(combined_hash, combined_hash, hash_Yi); //REF-40
    mpz_add(combined_hash, combined_hash, U_val->upki); //REF-40

    mpz_t n, d, e; //REF-01
    mpz_inits(n, d, e, NULL); //REF-02
    generate_rsa_keys(n, d, e, state);

    mpz_t signature; //REF-01
    mpz_init(signature); //REF-02
    mpz_init(U_val->signature); //REF-02
    mpz_init(rval->signature); //REF-02
    mpz_powm(signature, combined_hash, d, n); //REF-45

    mpz_set(U_val->rsa1, d); //REF-13
    mpz_set(U_val->rsa2, n); //REF-13
    mpz_set(U_val->signature, signature); //REF-13
    mpz_set(rval->signature, signature); //REF-13

    mpz_t verified_message; //REF-01
    mpz_init(verified_message); //REF-02
    mpz_powm(verified_message, signature, e, n); //REF-45
    mpz_set(rval->rsa2, n); //REF-13
    mpz_set(rval->rsa3, e); //REF-13

    printf("Group Manager GM processing : ");
    if (mpz_cmp(verified_message, combined_hash) == 0) //REF-07
        printf("\n----GM verifies U%lu digital Signature----\n", j);
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");

    element_t Xi_1; //REF-19
    element_init_G1(Xi_1, retval->pairing); //REF-20
    element_pow_zn(Xi_1, g1, xi); //REF-34
    element_printf("GM sends X%lu_1 = %B\n", j, Xi_1); //REF-30

    printf("User U%lu processing : ", j);
    element_pairing(temp1, Xi_1, retval->G2); //REF-38
    element_pairing(temp2, retval->G1, Xi_2); //REF-38
    if (element_cmp(temp1, temp2) == 0) //REF-45
        printf("\n----Verification passed : U%lu owns the membership----\n", j);
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");

    printf("User%lu Secret Key : \n", j);
    printf("   Valid Membership Certificates : \n");
    element_printf("\tA%lu = %B\n", j, Ai); //REF-30
    element_printf("\tX%lu = %B %B\n", j, Xi_1, Xi_2); //REF-30
    element_printf("\ty%lu = %B\n", j, yi); //REF-30
    printf("   Attribute certificate :\n");
    for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
        long unsigned int x = A[i];
        element_printf("  T%lu_%lu : %B\n", j, x, Ti[x]); //REF-30
    }

    element_init_Zr(U_val->yi, retval->pairing); //REF-31
    element_set(U_val->yi, yi); //REF-28
    element_init_G1(U_val->Ai, retval->pairing); //REF-20
    element_set(U_val->Ai, Ai); //REF-28
    element_init_G2(U_val->Xi_1, retval->pairing); //REF-21
    element_init_G2(U_val->Xi_2, retval->pairing); //REF-21
    element_set(U_val->Xi_1, Xi_1); //REF-28
    element_set(U_val->Xi_2, Xi_2); //REF-28
    element_init_G1(U_val->Yi, retval->pairing); //REF-20
    element_set(U_val->Yi, Yi); //REF-28
    element_init_Zr(U_val->xi, retval->pairing); //REF-31
    U_val->Ti = Ti;
    U_val->a = mpz_get_ui(a); //REF-09
    U_val->A = A;
    element_set(U_val->xi, xi); //REF-28

    element_init_G1(rval->Ai, retval->pairing); //REF-20
    element_set(rval->Ai, Ai); //REF-28
    element_init_G2(rval->Xi_1, retval->pairing); //REF-21
    element_init_G2(rval->Xi_2, retval->pairing); //REF-21
    element_set(rval->Xi_1, Xi_1); //REF-28
    element_set(rval->Xi_2, Xi_2); //REF-28
    element_init_G1(rval->Yi, retval->pairing); //REF-20
    element_set(rval->Yi, Yi); //REF-28

    printf("\n");
}

void buildtree(public *retval, mpz_t *S) {
    printf("\n----------------Build Tree Algorithm------------------\n");
    long unsigned int i = 0;
    mpz_t max, dummy, temp1; //REF-01
    mpz_inits(max, dummy, temp1, NULL); //REF-02
    while (1) {
        mpz_ui_pow_ui(temp1, 2, i); //REF-47
        if (mpz_cmp(retval->a, temp1) > 0) {
            //REF-11
            i++;
        } else
            break;
    }
    mpz_ui_pow_ui(dummy, 2, i); //REF-47
    mpz_sub_ui(dummy, dummy, 1); //REF-43
    mpz_ui_pow_ui(temp1, 2, i + 1); //REF-47
    mpz_sub_ui(temp1, temp1, 1); //REF-43
    mpz_add(max, temp1, dummy); //REF-40

    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node)); //REF-09
    Node *T_d = (Node *) malloc(mpz_get_ui(dummy) * sizeof(Node)); //REF-09
    for (i = 0; i < mpz_get_ui(temp1); i++) {
        //REF-09
        initializeNode(&T[i]);
    }
    for (i = 0; i < mpz_get_ui(dummy); i++) {
        //REF-09
        initializeNode(&T_d[i]);
    }

    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        //REF-09
        i = j - mpz_get_ui(dummy); //REF-09
        mpf_set_z(T[j].secret_value, S[i]); //REF-48
    }

    mpz_set_ui(T[0].name, 1); //REF-03
    long unsigned int k = 1, l = 0, m = 2;
    while (m <= mpz_get_ui(max)) {
        //REF-09
        mpz_set_ui(T[k].name, m); //REF-03
        mpz_set_ui(T[k + 1].name, m + 1); //REF-03
        mpz_set_ui(T_d[l].name, m + 2); //REF-03
        k = k + 2;
        l++;
        m = m + 3;
    }

    i = mpz_get_ui(temp1) - 1; //REF-09
    while (i > 0) {
        long unsigned int j = (i / 2) - 1;
        mpf_t temp4, temp5; //REF-49
        mpf_inits(temp4, temp5, NULL); //REF-50
        mpf_set_z(temp4, T[i - 1].name); //REF-48
        mpf_set_z(temp5, T[i].name); //REF-48
        mpf_mul(temp4, T[i].secret_value, temp4); //REF-51
        mpf_mul(temp5, T[i - 1].secret_value, temp5); //REF-51
        mpf_sub(T[j].secret_value, temp5, temp4); //REF-52

        mpf_mul_ui(temp4, T[i].secret_value, 2); //REF-53
        mpf_sub(T_d[j].secret_value, temp4, T[i - 1].secret_value); //REF-52
        i = i - 2;
    }
    l = 0, m = 1, k = 1;
    gmp_printf("Root  : %.Ff\n", T[0].secret_value); //REF-14
    while (m < mpz_get_ui(max)) {
        gmp_printf("Node %Zd: %.Ff\n", T[k].name, T[k].secret_value); //REF-14
        gmp_printf("Node %Zd: %.Ff\n", T[k + 1].name, T[k + 1].secret_value); //REF-14
        k = k + 2;
        gmp_printf("Node %Zd: %.Ff\n", T_d[l].name, T_d[l].secret_value); //REF-14
        l++;
        m = m + 3;
    }
    mpz_init_set(retval->max, max); //REF-54
    mpz_init_set(retval->temp1, temp1); //REF-54
    mpz_init_set(retval->d, dummy); //REF-54
    retval->T_d = T_d;
    private.T = T;
    element_init_G2(retval->vT, retval->pairing); //REF-21
    mpz_set_f(temp1, T[0].secret_value); //REF-55
    mpz_abs(temp1, temp1); //REF-56
    element_t tem1; //REF-19
    element_init_Zr(tem1, ret_setup.pairing); //REF-31
    element_set_mpz(tem1, temp1); //REF-57
    element_pow_zn(retval->vT, retval->G2, tem1); //REF-34
    printf("\n");
}

void buildtreevalidity(public *retval) {
    printf("\n----------------Build Tree Validity Algorithm------------------\n");
    gmp_randstate_t state; //REF-04
    gmp_randinit_default(state); //REF-05
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1)); //REF-06
    mpz_t a; //REF-01
    mpz_init(a); //REF-02
    while (mpz_get_ui(a) == 0) {
        //REF-09
        mpz_urandomm(a, state, retval->a); //REF-33
    }

    long unsigned int *A = (long unsigned int *) malloc(sizeof(long unsigned int) * mpz_get_ui(a));
    long unsigned int k = 0;
    printf("Attributes :");
    do {
        int z = 0;
        unsigned long int y = mpz_get_ui(retval->a);
        unsigned long int temp1 = rand() % y;
        for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
            if (temp1 == A[i]) {
                z = 1;
                break;
            }
        }
        if (z == 0) {
            A[k] = temp1;
            k++;
        }
    } while (k < mpz_get_ui(a)); //REF-09
    insertionSort(A, mpz_get_ui(a));
    for (k = 0; k < mpz_get_ui(a); k++) {
        gmp_printf("\t%lu", A[k] + 1); //REF-14
    }
    verify_tree(A, retval->d, retval->temp1, retval->max, ret_setup.T_d, 0);
}

void sign(user_output *U_val) {
    gmp_printf("\n*----------------Sign Algorithm for U%lu------------------\n", U_val->i); //REF-14
    printf("Attributes :");
    for (long unsigned int i = 0; i < U_val->a; i++) {
        gmp_printf("\t%lu", U_val->A[i] + 1); //REF-14
    }
    verify_tree(U_val->A, ret_setup.d, ret_setup.temp1, ret_setup.max, ret_setup.T_d, U_val->i);

    element_t temp, ID, z, r; //REF-19
    element_init_G1(ID, ret_setup.pairing); //REF-20
    element_init_Zr(temp, ret_setup.pairing); //REF-31
    element_init_Zr(z, ret_setup.pairing); //REF-31
    element_init_Zr(r, ret_setup.pairing); //REF-31
    element_init_Zr(R[U_val->i - 1].z, ret_setup.pairing); //REF-31
    element_init_Zr(R[U_val->i - 1].r, ret_setup.pairing); //REF-31

    element_random(z); //REF-25
    element_add(temp, z, U_val->yi); //REF-35
    element_invert(temp, temp); //REF-36
    element_pow_zn(ID, ret_setup.G1, temp); //REF-34
    element_random(r); //REF-25
    element_set(R[U_val->i - 1].r, r); //REF-23
    element_set(R[U_val->i - 1].z, z); //REF-23

    element_t *rho3 = (element_t *) malloc(sizeof(element_t) * 2); //REF-01
    element_t *rho7 = (element_t *) malloc(sizeof(element_t) * 2); //REF-01
    element_t rho4, rho5, rho6, rho8, rho9, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, vT; //REF-01
    element_init_G1(rho3[0], ret_setup.pairing); //REF-20
    element_init_G2(rho3[1], ret_setup.pairing); //REF-21
    element_init_G1(rho4, ret_setup.pairing); //REF-20
    element_init_G1(rho5, ret_setup.pairing); //REF-20
    element_init_G1(rho6, ret_setup.pairing); //REF-20
    element_init_G1(rho7[0], ret_setup.pairing); //REF-20
    element_init_G2(rho7[1], ret_setup.pairing); //REF-21
    element_init_G1(rho8, ret_setup.pairing); //REF-20
    element_init_G2(rho9, ret_setup.pairing); //REF-21
    element_init_Zr(temp1, ret_setup.pairing); //REF-31
    element_init_Zr(temp2, ret_setup.pairing); //REF-31
    element_init_G1(temp3, ret_setup.pairing); //REF-20
    element_init_G1(temp4, ret_setup.pairing); //REF-20
    element_init_G2(temp5, ret_setup.pairing); //REF-21
    element_init_G2(temp9, ret_setup.pairing); //REF-21
    element_init_GT(temp6, ret_setup.pairing); //REF-22
    element_init_GT(temp7, ret_setup.pairing); //REF-22
    element_init_GT(temp8, ret_setup.pairing); //REF-22
    element_init_G2(vT, ret_setup.pairing); //REF-21
    element_init_G2(U_val->vT, ret_setup.pairing); //REF-21

    element_pow_mpz(vT, ret_setup.G2, R[U_val->i - 1].sT); //REF-26
    element_set(U_val->vT, vT); //REF-23

    element_set(rho3[0], U_val->Xi_1); //REF-23
    element_set(rho3[1], U_val->Xi_2); //REF-23

    element_add(temp1, private.gamma, U_val->xi); //REF-35
    element_invert(temp2, temp1); //REF-36
    if (mpz_cmp_ui(R[U_val->i - 1].sT1, 0) > 0) {
        //REF-07
        element_mul_mpz(temp1, temp2, R[U_val->i - 1].sT1); //REF-37
        element_pow_zn(rho4, ret_setup.generators_G1[0], temp1); //REF-34
    } else {
        mpz_t tem1; //REF-01
        mpz_init(tem1); //REF-23
        mpz_abs(tem1, R[U_val->i - 1].sT1); //REF-56
        element_mul_mpz(temp2, temp2, tem1); //REF-37
        element_pow_zn(rho4, ret_setup.generators_G1[0], temp2); //REF-34
        element_invert(rho4, rho4); //REF-36
    }

    if (mpz_cmp_ui(R[U_val->i - 1].sT2, 0) > 0) {
        //REF-07
        element_pow_mpz(rho5, ret_setup.generators_G1[0], R[U_val->i - 1].sT2); //REF-26
    } else {
        mpz_t tem1; //REF-01
        mpz_init_set(tem1, R[U_val->i - 1].sT2); //REF-54
        mpz_abs(tem1, tem1); //REF-56
        element_pow_mpz(rho5, ret_setup.generators_G1[0], tem1); //REF-26
        element_invert(rho5, rho5); //REF-36
    }

    element_set(rho6, ID); //REF-23

    element_pow_zn(rho7[0], ret_setup.G1, z); //REF-34
    element_pow_zn(rho7[1], ret_setup.G2, z); //REF-34

    element_set(temp3, ret_setup.wf); //REF-23
    element_pow_zn(temp3, temp3, r); //REF-34
    element_pow_zn(temp4, ret_setup.generators_G1[0], z); //REF-34
    element_mul(rho8, temp3, temp4); //REF-29

    element_pow_zn(rho9, ret_setup.G2, r); //REF-34

    element_mul(temp5, ret_setup.omega, rho3[1]); //REF-29
    element_pairing(temp6, U_val->Ai, temp5); //REF-38
    element_pairing(temp7, ret_setup.generators_G1[0], ret_setup.G2); //REF-38
    element_pow_zn(temp9, ret_setup.G2, U_val->yi); //REF-34
    element_pairing(temp8, ret_setup.G1, temp9); //REF-38
    element_mul(temp7, temp7, temp8); //REF-29
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----Signer has valid Membership Cerficate-----");
    else
        printf("-*-*-Fraud-*-*-");

    element_pairing(temp6, rho4, temp5); //REF-38
    element_pairing(temp7, rho5, ret_setup.G2); //REF-38
    element_pairing(temp8, ret_setup.generators_G1[0], U_val->vT); //REF-38
    element_mul(temp6, temp6, temp7); //REF-29
    if (element_cmp(temp6, temp8) == 0) //REF-27
        printf("\n-----Signer has valid Attribute Cerficate-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pow_zn(temp5, ret_setup.G2, U_val->yi); //REF-34
    element_mul(temp9, temp5, rho7[1]); //REF-29
    element_pairing(temp6, rho6, temp9); //REF-38
    element_pairing(temp7, ret_setup.G1, ret_setup.G2); //REF-38
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----Signer has valid ID-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho8, ret_setup.G2); //REF-38
    element_pairing(temp7, ret_setup.generators_G1[0], rho7[1]); //REF-38
    element_pairing(temp8, ret_setup.wf, rho9); //REF-38
    element_mul(temp7, temp7, temp8); //REF-29
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----Water Signature valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pow_zn(temp4, ret_setup.G1, U_val->yi); //REF-34
    element_pairing(temp6, temp4, ret_setup.G2); //REF-38
    element_pow_zn(temp5, ret_setup.G2, U_val->yi); //REF-34
    element_pairing(temp7, ret_setup.G1, temp5); //REF-38
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----yi commited valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho3[0], ret_setup.G2); //REF-38
    element_pairing(temp7, ret_setup.G1, rho3[1]); //REF-38
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----Xi commited valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho7[0], ret_setup.G2); //REF-38
    element_pairing(temp7, ret_setup.G1, rho7[1]); //REF-38
    if (element_cmp(temp6, temp7) == 0) //REF-27
        printf("\n-----Non-frameability Adversary checked valid-----\n");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_init_G1(U_val->rho1, ret_setup.pairing); //REF-20
    element_init_Zr(U_val->rho2, ret_setup.pairing); //REF-31
    element_init_G1(U_val->rho4, ret_setup.pairing); //REF-20
    element_init_G1(U_val->rho5, ret_setup.pairing); //REF-20
    element_init_G1(U_val->rho6, ret_setup.pairing); //REF-20
    element_init_G1(U_val->rho8, ret_setup.pairing); //REF-20
    element_init_G2(U_val->rho9, ret_setup.pairing); //REF-21
    element_set(U_val->rho1, U_val->Ai); //REF-28
    element_set(U_val->rho2, U_val->yi); //REF-28
    element_set(U_val->rho4, rho4); //REF-28
    element_set(U_val->rho5, rho5); //REF-28
    element_set(U_val->rho6, rho6); //REF-28
    element_set(U_val->rho8, rho8); //REF-28
    element_set(U_val->rho9, rho9); //REF-28
    U_val->rho3 = rho3;
    U_val->rho7 = rho7;

    element_printf("\n  User U%lu Signature :", U_val->i); //REF-30
    element_printf("\n\trho1 : %B", U_val->Ai); //REF-30
    element_printf("\n\trho2 : %B", U_val->yi); //REF-30
    element_printf("\n\trho3,1 : %B", rho3[0]); //REF-30
    element_printf("\n\trho3,2 : %B", rho3[1]); //REF-30
    element_printf("\n\trho4 : %B", rho4); //REF-30
    element_printf("\n\trho5 : %B", rho5); //REF-30
    element_printf("\n\trho6 : %B", rho6); //REF-30
    element_printf("\n\trho7,1 : %B", rho7[0]); //REF-30
    element_printf("\n\trho7,2 : %B", rho7[1]); //REF-30
    element_printf("\n\trho8 : %B", rho8); //REF-30
    element_printf("\n\trho9 : %B", rho9); //REF-30

    element_t *sigma1 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma2_1 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma2_2 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma3_1 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma3_2 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma4 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
    element_t *sigma5 = (element_t *) malloc(sizeof(element_t) * 2); //REF-19

    element_t tem1[2], tem2[2], tem3[2], tem4[2], tem5[2], tem6[2], tem7[2], tem8[2]; //REF-19
    for (int i = 0; i < 2; i++) {
        element_init_G1(tem1[i], ret_setup.pairing); //REF-20
        element_init_G1(tem2[i], ret_setup.pairing); //REF-20
        element_init_G1(tem3[i], ret_setup.pairing); //REF-20
        element_init_G1(tem4[i], ret_setup.pairing); //REF-20
        element_init_G2(tem5[i], ret_setup.pairing); //REF-21
        element_init_G2(tem6[i], ret_setup.pairing); //REF-21
        element_init_G2(tem7[i], ret_setup.pairing); //REF-21
        element_init_G1(sigma1[i], ret_setup.pairing); //REF-20
        element_init_G1(sigma2_1[i], ret_setup.pairing); //REF-21
        element_init_G2(sigma2_2[i], ret_setup.pairing); //REF-21
        element_init_G1(sigma3_1[i], ret_setup.pairing); //REF-20
        element_init_G2(sigma3_2[i], ret_setup.pairing); //REF-21
        element_init_G1(sigma4[i], ret_setup.pairing); //REF-20
        element_init_G1(sigma5[i], ret_setup.pairing); //REF-20
    }

    element_init_G1(tem8[0], ret_setup.pairing); //REF-20
    element_init_G2(tem8[1], ret_setup.pairing); //REF-21
    element_set_si(tem8[0], 1); //REF-58
    element_set_si(tem8[1], 1); //REF-58

    element_pow_zn(tem1[0], ret_setup.u[0], r); //REF-34
    element_pow_zn(tem1[1], ret_setup.u[1], r); //REF-34
    element_pow_zn(tem2[0], ret_setup.u[2], z); //REF-34
    element_pow_zn(tem2[1], ret_setup.u[3], z); //REF-34
    element_mul(tem4[0], tem2[0], tem1[0]); //REF-29
    element_mul(tem4[1], tem2[1], tem1[1]); //REF-29
    element_mul(sigma1[0], tem4[0], tem8[0]); //REF-29
    element_mul(sigma1[1], tem4[1], U_val->rho1); //REF-29
    element_printf("\n\tsigma1 : %B %B", sigma1[0], sigma1[1]); //REF-30

    element_mul(tem2[0], ret_setup.u[2], tem8[0]); //REF-29
    element_mul(tem2[1], ret_setup.u[3], ret_setup.G1); //REF-29
    element_pow_zn(tem3[0], tem2[0], U_val->rho2); //REF-34
    element_pow_zn(tem3[1], tem2[1], U_val->rho2); //REF-34
    element_mul(sigma2_1[0], tem3[0], tem1[0]); //REF-29
    element_mul(sigma2_1[1], tem3[1], tem2[1]); //REF-29
    element_printf("\n\tsigma2_1 : %B %B", sigma2_1[0], sigma2_1[1]); //REF-30

    element_pow_zn(tem5[0], ret_setup.v[0], r); //REF-34
    element_pow_zn(tem5[1], ret_setup.v[1], r); //REF-34
    element_mul(tem6[0], ret_setup.v[2], tem8[1]); //REF-29
    element_mul(tem6[1], ret_setup.v[3], ret_setup.G2); //REF-29
    element_pow_zn(tem7[0], tem6[0], U_val->rho2); //REF-34
    element_pow_zn(tem7[1], tem6[1], U_val->rho2); //REF-34
    element_mul(sigma2_2[0], tem7[0], tem5[0]); //REF-29
    element_mul(sigma2_2[1], tem7[1], tem5[1]); //REF-29
    element_printf("\n\tsigma2_2 : %B %B", sigma2_2[0], sigma2_2[1]); //REF-30

    element_mul(sigma3_1[0], tem4[0], tem8[0]); //REF-29
    element_mul(sigma3_1[1], tem4[1], U_val->rho3[0]); //REF-29
    element_printf("\n\tsigma3,1 : %B %B", sigma3_1[0], sigma3_1[1]); //REF-30
    element_pow_zn(tem5[0], ret_setup.v[0], r); //REF-34
    element_pow_zn(tem5[1], ret_setup.v[1], r); //REF-34
    element_pow_zn(tem6[0], ret_setup.v[2], z); //REF-34
    element_pow_zn(tem6[1], ret_setup.v[3], z); //REF-34
    element_mul(tem6[0], tem6[0], tem5[0]); //REF-29
    element_mul(tem6[1], tem6[1], tem5[1]); //REF-29
    element_mul(sigma3_2[0], tem6[0], tem8[1]); //REF-29
    element_mul(sigma3_2[1], tem6[1], U_val->rho3[1]); //REF-29
    element_printf("\n\tsigma3,2 : %B %B", sigma3_2[0], sigma3_2[1]); //REF-30

    element_pow_zn(tem1[0], ret_setup.u_dash[0], r); //REF-34
    element_pow_zn(tem1[1], ret_setup.u_dash[1], r); //REF-34
    element_pow_zn(tem2[0], ret_setup.u_dash[2], z); //REF-34
    element_pow_zn(tem2[1], ret_setup.u_dash[3], z); //REF-34
    element_mul(tem4[0], tem2[0], tem1[0]); //REF-29
    element_mul(tem4[1], tem2[1], tem1[1]); //REF-29
    element_mul(sigma4[0], tem4[0], tem8[0]); //REF-29
    element_mul(sigma4[1], tem4[1], U_val->rho4); //REF-29
    element_printf("\n\tsigma4 : %B %B", sigma4[0], sigma4[1]); //REF-30

    element_mul(sigma5[0], tem4[0], tem8[0]); //REF-29
    element_mul(sigma5[1], tem4[1], U_val->rho5); //REF-29
    element_printf("\n\tsigma5 : %B %B\n", sigma5[0], sigma5[1]); //REF-30

    U_val->sigma1 = sigma1;
    U_val->sigma2_1 = sigma2_1;
    U_val->sigma2_2 = sigma2_2;
    U_val->sigma3_1 = sigma3_1;
    U_val->sigma3_2 = sigma3_2;
    U_val->sigma4 = sigma4;
    U_val->sigma5 = sigma5;
}

void verify(user_output *U_val) {
    mpz_t temp1; //REF-01
    mpz_init(temp1); //REF-02
    element_t groot, vT, tem1, tem2, sT; //REF-19
    element_init_G1(groot, ret_setup.pairing); //REF-20
    element_init_G2(vT, ret_setup.pairing); //REF-21
    element_init_GT(tem1, ret_setup.pairing); //REF-22
    element_init_GT(tem2, ret_setup.pairing); //REF-20
    element_init_Zr(sT, ret_setup.pairing); //REF-31
    mpz_set_f(temp1, U_val->T[0].secret_value); //REF-55
    mpz_abs(temp1, temp1); //REF-56
    element_set_mpz(sT, temp1); //REF-57
    element_pow_zn(groot, ret_setup.G1, sT); //REF-34
    element_set(vT, ret_setup.vT); //REF-28
    element_pairing(tem1, groot, ret_setup.G2); //REF-38
    element_pairing(tem2, ret_setup.G1, vT); //REF-38
    if (element_cmp(tem1, tem2) == 0) //REF-27
        gmp_printf("\n-----Verify algorithm valid for User%lu-----\n", U_val->i); //REF-14

    else
        printf("\nx-x-x-Verify algorithm invalid-x-x-x\n");
}

void openuser(user_output *U_val, reg *rval, element_t ok_user) {
    if (private.alpha1 == ok_user) {
        printf("\n*----------OpenUser verified-----------");

        mpz_t hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki; //REF-01
        mpz_inits(hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki, NULL); //REF-02
        element_to_mpz(hash_Ai, U_val->Ai); //REF-39
        element_to_mpz(hash_Xi_2, U_val->Xi_2); //REF-39
        element_to_mpz(hash_Yi, U_val->Yi); //REF-39

        mpz_add(combined_hash, hash_Ai, hash_Xi_2); //REF-40
        mpz_add(combined_hash, combined_hash, hash_Yi); //REF-40
        mpz_add(combined_hash, combined_hash, rval->upki); //REF-40
        mpz_t verified_message; //REF-01
        mpz_init(verified_message); //REF-02
        mpz_powm(verified_message, rval->signature, rval->rsa3, rval->rsa2); //REF-45
        if (mpz_cmp(verified_message, combined_hash) == 0) {
            //REF-11
            printf("\n  ----OpenUser verifies U%lu digital Signature----\n", rval->i);

            element_t *sigma1_verify = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
            element_t tem1[2], tem2[2], tem4[2]; //REF-19
            for (int i = 0; i < 2; i++) {
                element_init_G1(tem1[i], ret_setup.pairing); //REF-20
                element_init_G1(tem2[i], ret_setup.pairing); //REF-20
                element_init_G1(tem4[i], ret_setup.pairing); //REF-20
                element_init_G1(sigma1_verify[i], ret_setup.pairing); //REF-20
            }
            element_pow_zn(tem1[0], ret_setup.u[0], rval->r); //REF-34
            element_pow_zn(tem1[1], ret_setup.u[1], rval->r); //REF-34
            element_pow_zn(tem2[0], ret_setup.u[2], rval->z); //REF-34
            element_pow_zn(tem2[1], ret_setup.u[3], rval->z); //REF-34
            element_mul(tem4[0], tem2[0], tem1[0]); //REF-29
            element_mul(tem4[1], tem2[1], tem1[1]); //REF-29
            element_div(sigma1_verify[0], U_val->sigma1[0], tem4[0]); //REF-46
            element_div(sigma1_verify[1], U_val->sigma1[1], tem4[1]); //REF-46

            if (element_cmp(rval->Ai, sigma1_verify[1]) == 0) {
                //REF-27
                gmp_printf("\t---User Identity %lu ---\n", rval->i); //REF-14
            } else
                printf("\n\t-*-*-Wrong or No User Identity-*-*- ");
        } else
            printf("\n-x-x-x-Verification failed-x-x-x-\n");
    } else
        printf("\n-x-x-x-Open User Verification failed-x-x-x-\n");
}

void traceatt(user_output *U_val, reg *rval, element_t tk_att) {
    if (private.alpha1_dash == tk_att) {
        gmp_printf("\n----Tracer is valid----\n"); //REF-14

        element_t *sigma5_verify = (element_t *) malloc(sizeof(element_t) * 2); //REF-19
        element_t tem1[2], tem2[2], tem4[2]; //REF-19
        for (int i = 0; i < 2; i++) {
            element_init_G1(tem1[i], ret_setup.pairing); //REF-20
            element_init_G1(tem2[i], ret_setup.pairing); //REF-20
            element_init_G1(tem4[i], ret_setup.pairing); //REF-20
            element_init_G1(sigma5_verify[i], ret_setup.pairing); //REF-20
        }
        element_pow_zn(tem1[0], ret_setup.u_dash[0], rval->r); //REF-34
        element_pow_zn(tem1[1], ret_setup.u_dash[1], rval->r); //REF-34
        element_pow_zn(tem2[0], ret_setup.u_dash[2], rval->z); //REF-34
        element_pow_zn(tem2[1], ret_setup.u_dash[3], rval->z); //REF-34
        element_mul(tem4[0], tem2[0], tem1[0]); //REF-29
        element_mul(tem4[1], tem2[1], tem1[1]); //REF-29
        element_div(sigma5_verify[0], U_val->sigma5[0], tem4[0]); //REF-46
        element_div(sigma5_verify[1], U_val->sigma5[1], tem4[1]); //REF-46

        element_t rho5_verify; //REF-19
        element_init_G1(rho5_verify, ret_setup.pairing); //REF-20
        if (mpz_cmp_ui(rval[U_val->i - 1].sT2, 0) > 0) {
            //REF-07
            element_pow_mpz(rho5_verify, ret_setup.generators_G1[0], R[U_val->i - 1].sT2); //REF-26
        } else {
            mpz_t temp1; //REF-01
            mpz_init_set(temp1, R[U_val->i - 1].sT2); //REF-54
            mpz_abs(temp1, temp1); //REF-56
            element_pow_mpz(rho5_verify, ret_setup.generators_G1[0], temp1); //REF-26
            element_invert(rho5_verify, rho5_verify); //REF-36
        }

        gmp_printf("U%lu Attribute Name: ", U_val->i);
        if (element_cmp(rho5_verify, sigma5_verify[1]) == 0) {
            //REF-27
            for (long unsigned int i = 0; i < (U_val->a); i++) {
                gmp_printf("   %lu", U_val->A[i] + 1); //REF-14
            }
            printf("\n");
        } else
            printf("\n\t-*-*-Wrong or No User Identity-*-*- ");
    } else
        gmp_printf("\n----Tracer is invalid----\n"); //REF-14
}

int main() {
    srand(time(NULL)); // define the state for random functions
    mpz_t security_parameter; //REF-01
    mpz_init(security_parameter); //REF-02
    mpz_set_ui(security_parameter, 16); //REF-03

    setup(&ret_setup, security_parameter);

    // Generate a random message of O(security_parameter) bits
    long unsigned int m_dash = mpz_get_ui(security_parameter);
    char *message = (char *) malloc(sizeof(char) * m_dash);
    for (unsigned long int i = 0; i < m_dash; i++) {
        if ((rand() % 2) == 0)
            message[i] = '0';
        else
            message[i] = '1';
    }
    message[m_dash] = '\0';
    ret_setup.message = message;
    printf("Message (m' bits): %s\n", message);
    ret_setup.message = message;
    keygen(message, &ret_setup, &private);

    for (int i = 0; i < num_user; i++) {
        join(&ret_setup, i + 1, &U[i], &R[i], &private);
        insertionSort(U[i].A, U[i].a);
    }

    buildtree(&ret_setup, private.S);

    buildtreevalidity(&ret_setup);

    for (int i = 0; i < num_user; i++) {
        sign(&U[i]);
        verify(&U[i]);
    }

    for (int i = 0; i < num_user; i++) {
        openuser(&U[i], &R[i], private.alpha1);
    }

    for (int i = 0; i < num_user; i++) {
        traceatt(&U[i], &R[i], private.alpha1_dash);
    }

    printf("\n'::+-'*-++-*'-+::'");
    return 0;
}
