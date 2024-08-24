#include <stdio.h>
#include <pbc/pbc.h>
#include <gmp.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define num_user 1

typedef struct Node {
    mpz_t name;
    mpf_t secret_value;
    long unsigned int threshold;
} Node;

typedef struct public {
    pbc_param_t param;
    element_t G1;
    element_t G2;
    mpz_t p;
    mpz_t k;
    mpz_t a;
    pairing_t pairing;
    element_t identity_g1;
    element_t identity_g2;
    element_t *generators_G1;
    element_t *generators_G2;

    char *message;
    element_t wf;
    element_t *gatt;
    element_t omega;
    element_t *u;
    element_t *u_dash;
    element_t *v;

    element_t vT;
    mpz_t d;
    Node *T_d;
    Node *T_trial;
    mpz_t max;
    mpz_t temp1;
    mpz_t dummy;
}
public;

public
ret_setup;

typedef struct secret {
    mpz_t *S;
    element_t gamma;
    element_t alpha1;
    element_t alpha1_dash;
    Node *T;
} secret;

secret private;

typedef struct user_output {
    long unsigned int i;
    mpz_t upki;
    mpz_t uski;
    element_t Ai;
    element_t Xi_1;
    element_t Xi_2;
    element_t yi_dash_dash;

    mpz_t signature;
    mpz_t rsa1;
    mpz_t rsa2;
    element_t Yi;
    element_t yi;
    element_t *Ti;
    element_t xi;
    long unsigned int a;
    long unsigned int *A;
    mpz_t d;
    long unsigned int *D;

    Node *T;
    Node *T_d;

    element_t vT;

    element_t rho1;
    element_t rho2;
    element_t *rho3;
    element_t rho4;
    element_t rho5;
    element_t rho6;
    element_t *rho7;
    element_t rho8;
    element_t rho9;
    element_t *sigma1;
    element_t *sigma2_1;
    element_t *sigma2_2;
    element_t *sigma3_1;
    element_t *sigma3_2;
    element_t *sigma4;
    element_t *sigma5;
} user_output;

user_output U[num_user];

typedef struct reg {
    long unsigned int i;
    mpz_t upki;
    long unsigned int a;
    long unsigned int *A;
    mpz_t d;
    long unsigned int *D;
    element_t Ai;
    element_t Xi_1;
    element_t Xi_2;
    element_t Yi;
    mpz_t signature;
    mpz_t rsa2;
    mpz_t rsa3;

    mpz_t sT;
    mpz_t sT1;
    mpz_t sT2;
    element_t r;
    element_t z;
} reg;

reg R[num_user];

void random_prime_bits(mpz_t result, mpz_t n) {
    mpz_t random;
    mpz_init(random);

    gmp_randstate_t state;

    gmp_randinit_default(state);

    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1));

    if (mpz_cmp_ui(n, 1) <= 0) {
        printf("NO PRIME EXISTS\n");
    } else {
        mpz_t lower_limit;
        mpz_init(lower_limit);
        mpz_ui_pow_ui(lower_limit, 2, mpz_get_ui(n) - 1);

        while (1) {
            mpz_urandomb(random, state, mpz_get_ui(n));

            if (mpz_cmp(random, lower_limit) > 0 && mpz_probab_prime_p(random, mpz_get_ui(n))) {
                mpz_set(result, random);
                break;
            }
        }
    }
}

void generate_rsa_keys(mpz_t n, mpz_t d, mpz_t e, gmp_randstate_t state) {
    mpz_t p, q, phi;
    mpz_inits(p, q, phi, NULL);
    mpz_set(e, ret_setup.p);

    mpz_urandomb(p, state, mpz_get_ui(ret_setup.k));
    mpz_urandomb(q, state, mpz_get_ui(ret_setup.k));
    mpz_nextprime(p, p);
    mpz_nextprime(q, q);

    mpz_mul(n, p, q);

    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_mul(phi, p, q);

    mpz_invert(d, e, phi);

    mpz_clears(p, q, phi, NULL);
}

void insertionSort(long unsigned int arr[], long unsigned int n) {
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void initializeNode(Node *node) {
    mpz_init(node->name);
    mpf_init(node->secret_value);
    node->threshold = 1;
}

void tree_making(mpz_t dummy, mpz_t temp1, mpz_t max, long unsigned int A[], Node *T_d, mpz_t *S, long unsigned int z) {
    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node));
    long unsigned int i = 0;
    for (i = 0; i < mpz_get_ui(temp1); i++) {
        initializeNode(&T[i]);
    }
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        i = j - mpz_get_ui(dummy);
        mpf_set_z(T[j].secret_value, S[i]);
    }
    mpz_set_ui(T[0].name, 1);
    long unsigned int k = 1, m = 2, l = 0;
    while (m <= mpz_get_ui(max)) {
        mpz_set_ui(T[k].name, m);
        mpz_set_ui(T[k + 1].name, m + 1);
        mpz_set_ui(T_d[l].name, m + 2);
        k = k + 2;
        l++;
        m = m + 3;
    }
    i = 0;
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        if (j - mpz_get_ui(dummy) == A[i]) {
            i++;
        } else {
            initializeNode(&T[j]);
        }
    }
    i = mpz_get_ui(temp1) - 1;
    printf("\n");
    mpf_t temp2, temp3, temp4;
    mpf_inits(temp2, temp3, temp4, NULL);
    mpz_t temp5;
    mpz_init(temp5);
    while (i > 0) {
        if (mpz_cmp_ui(T[i].name, 0) == 0 && mpz_cmp_ui(T[i - 1].name, 0) == 0) {
            initializeNode(&T[i / 2 - 1]);
            initializeNode(&T_d[i / 2 - 1]);
        } else if (mpz_cmp_ui(T[i - 1].name, 0) == 0) {
            mpf_set_z(temp2, T_d[i / 2 - 1].name);
            mpf_mul(temp2, temp2, T[i].secret_value);
            mpf_set_z(temp3, T[i].name);
            mpf_mul(temp3, temp3, T_d[i / 2 - 1].secret_value);
            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3);
        } else if (mpz_cmp_ui(T[i].name, 0) == 0) {
            mpf_set_z(temp2, T_d[i / 2 - 1].name);
            mpf_div_ui(temp2, temp2, 2);
            mpf_set_z(temp3, T[i - 1].name);
            mpf_div_ui(temp3, temp3, 2);
            mpf_mul(temp2, T[i - 1].secret_value, temp2);
            mpf_mul(temp3, T_d[i / 2 - 1].secret_value, temp3);
            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3);
        } else {
            mpz_mul(temp5, T[i - 1].name, T[i].name);
            mpf_set_z(temp2, temp5);
            mpf_mul(temp2, T_d[i / 2 - 1].secret_value, temp2);

            mpz_mul(temp5, T[i].name, T_d[i / 2 - 1].name);
            mpf_set_z(temp3, temp5);
            mpf_mul(temp4, temp3, T[i - 1].secret_value);

            mpf_add(temp2, temp2, temp4);
            mpf_div_ui(temp2, temp2, 2);

            mpz_mul(temp5, T_d[i / 2 - 1].name, T[i - 1].name);
            mpf_set_z(temp3, temp5);
            mpf_mul(temp3, temp3, T[i].secret_value);

            mpf_sub(T[i / 2 - 1].secret_value, temp2, temp3);
        }
        i = i - 2;
    }
    if (z == 0)
        ret_setup.T_trial = T;
    else
        U[z - 1].T = T;

    i = mpz_get_ui(temp1) - 1;
    while (i > 0) {
        long unsigned int j = (i / 2) - 1;
        mpf_mul_ui(temp4, private.T[i].secret_value, 2);
        mpf_sub(T_d[j].secret_value, temp4, private.T[i - 1].secret_value);
        i = i - 2;
    }
    m = 2, l = 0;
    while (m <= mpz_get_ui(max)) {
        mpz_set_ui(T_d[l].name, m + 2);
        l++;
        m = m + 3;
    }
}

void verify_tree(long unsigned int A[], mpz_t dummy, mpz_t temp1, mpz_t max, Node *T_d, unsigned long int z) {
    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node));
    mpz_t *S = (mpz_t *) malloc(mpz_get_ui(ret_setup.a) * sizeof(mpz_t));
    long unsigned int i = 0;
    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        long unsigned int k = j - mpz_get_ui(dummy);
        if (k == A[i]) {
            i++;
            mpz_init_set(S[k], private.S[k]);
        } else
            mpz_init(S[k]);
    }
    tree_making(dummy, temp1, max, A, T_d, S, z);
    if (z == 0)
        T = ret_setup.T_trial;
    else
        T = U[z - 1].T;
    gmp_printf("Root  : %.Ff\n", T[0].secret_value);
    long unsigned int l = 0, m = 1, k = 1;
    while (m < mpz_get_ui(max)) {
        gmp_printf("Node %Zd: %.Ff\n", T[k].name, T[k].secret_value);
        gmp_printf("Node %Zd: %.Ff\n", T[k + 1].name, T[k + 1].secret_value);
        k = k + 2;
        gmp_printf("Node %Zd: %.Ff\n", T_d[l].name, T_d[l].secret_value);
        l++;
        m = m + 3;
    }
    mpz_t temp5;
    mpz_init(temp5);
    element_t groot, vT, tem1, tem2, sT;
    element_init_G1(groot, ret_setup.pairing);
    element_init_G2(vT, ret_setup.pairing);
    element_init_GT(tem1, ret_setup.pairing);
    element_init_GT(tem2, ret_setup.pairing);
    element_init_Zr(sT, ret_setup.pairing);
    mpz_set_f(temp5, T[0].secret_value);
    mpz_abs(temp5, temp5);
    element_set_mpz(sT, temp5);
    element_pow_zn(groot, ret_setup.G1, sT);
    element_set(vT, ret_setup.vT);
    element_pairing(tem1, groot, ret_setup.G2);
    element_pairing(tem2, ret_setup.G1, vT);
    if (element_cmp(tem1, tem2) == 0) {
        if (z != 0) {
            printf("\n-----Signing verification valid-----\n");
        } else
            printf("\n-----Buildtree algorithm valid-----\n");
    } else
        printf("\nx-x-x-Buildtree algorithm invalid-x-x-x\n");

    Node *T_dtrial = (Node *) malloc(mpz_get_ui(dummy) * sizeof(Node));
    for (i = 0; i < mpz_get_ui(dummy); i++) {
        initializeNode(&T_dtrial[i]);
    }
    tree_making(dummy, temp1, max, A, T_dtrial, S, 0);

    mpz_t sT1, sT2, s;
    mpz_inits(sT1, sT2, s, NULL);
    mpz_set_f(sT1, ret_setup.T_trial[0].secret_value);
    for (i = 0; i < mpz_get_ui(ret_setup.a); i++) {
        mpz_set_ui(S[i], 0);
    }

    tree_making(dummy, temp1, max, A, T_d, S, 0);
    mpz_set_f(sT2, ret_setup.T_trial[0].secret_value);
    mpz_add(s, sT1, sT2);
    if (z != 0) {
        mpz_inits(R[z - 1].sT, R[z - 1].sT1, R[z - 1].sT2, NULL);
        if (mpz_cmp_ui(s, 0) < 0) {
            mpz_abs(R[z - 1].sT, s);
            if (mpz_cmp_ui(sT1, 0) < 0) {
                mpz_abs(R[z - 1].sT1, sT1);
                mpz_neg(R[z - 1].sT2, sT2);
            } else {
                mpz_abs(R[z - 1].sT2, sT2);
                mpz_neg(R[z - 1].sT1, sT1);
            }
        } else {
            mpz_set(R[z - 1].sT1, sT1);
            mpz_set(R[z - 1].sT, s);
            mpz_set(R[z - 1].sT2, sT2);
        }
    }
    gmp_printf("s  : %Zd\n", s);
    gmp_printf("st1  : %Zd\n", sT1);
    gmp_printf("st2  : %Zd\n", sT2);
}

void setup(public *retval, mpz_t k) {
    printf("\n----------------Setup Algorithm------------------\n");
    mpz_t p;
    mpz_init(p);
    mpz_init(retval->a);
    mpz_set(retval->a, k);
    random_prime_bits(p, k);
    gmp_printf("P(order of bilinear group) = %Zd\n", p);
    gmp_printf("Total Attributes = %Zd\n", k);

    pairing_t pairing;
    pbc_param_t param;
    pbc_param_init_a1_gen(param, p);
    pairing_init_pbc_param(pairing, param);

    element_t g1, g2, gt1, identity_g1, identity_g2, temp_g1, temp_g2;
    element_init_G1(g1, pairing);
    element_init_G2(g2, pairing);
    element_init_GT(gt1, pairing);
    element_init_G1(temp_g1, pairing);
    element_init_G2(temp_g2, pairing);
    element_init_G1(identity_g1, pairing);
    element_init_G2(identity_g2, pairing);

    element_set0(identity_g2);
    element_set0(identity_g1);

    mpz_t required, gen;
    mpz_init(required);
    mpz_init(gen);
    mpz_add_ui(required, k, 2);

    element_t *generators_g1 = (element_t *) malloc(sizeof(element_t) * (mpz_get_ui(required)));
    element_t *generators_g2 = (element_t *) malloc(sizeof(element_t) * (mpz_get_ui(required)));

    for (unsigned long int i = 0; i < mpz_get_ui(required); i++) {
        element_init_G1(generators_g1[i], pairing);
        element_init_G2(generators_g2[i], pairing);
    }
    unsigned long long int index = 0;
    do {
        element_random(g1);
        element_pow_mpz(temp_g1, g1, p);
        if (element_cmp(temp_g1, identity_g1) == 0) {
            mpz_add_ui(gen, gen, 1);
            element_set(generators_g1[index], g1);
            index++;
        }
    } while (mpz_cmp(gen, required));
    index = 0;

    mpz_set_ui(gen, 0);
    do {
        element_random(g2);
        element_pow_mpz(temp_g2, g2, p);
        if (element_cmp(temp_g2, identity_g2) == 0) {
            mpz_add_ui(gen, gen, 1);
            element_set(generators_g2[index], g2);
            index++;
        }
    } while (mpz_cmp(gen, required));

    pbc_param_init_a1_gen(retval->param, p);
    pairing_init_pbc_param(retval->pairing, param);
    mpz_init(retval->k);
    mpz_set(retval->k, k);
    mpz_init(retval->p);
    mpz_set_ui(retval->p, mpz_get_ui(p));
    element_init_G1(retval->G1, pairing);
    element_init_G2(retval->G2, pairing);
    element_init_G1(retval->identity_g1, pairing);
    element_init_G2(retval->identity_g2, pairing);
    element_set0(retval->identity_g1);
    element_set0(retval->identity_g2);
    retval->generators_G1 = generators_g1;
    retval->generators_G2 = generators_g2;
}

void keygen(char *M, public *ret_val, secret *sec_val) {
    printf("\n----------------KenGen Algorithm------------------\n");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1));

    element_t wf;
    element_init_G1(wf, ret_val->pairing);
    element_set(wf, ret_val->generators_G1[1]);
    element_t temp_pw;
    element_init_G1(temp_pw, ret_val->pairing);
    mpz_t msg_val;
    mpz_init(msg_val);

    for (unsigned long int i = 0; i < mpz_get_ui(ret_val->k); i++) {
        if (M[i] == '1') {
            mpz_set_ui(msg_val, 1);
        } else {
            mpz_set_ui(msg_val, 0);
        }
        element_pow_mpz(temp_pw, ret_setup.generators_G1[i + 2], msg_val);
        element_mul(wf, wf, temp_pw);
    }

    printf("Group Public Key : \n");

    element_printf("  h = %B\n", ret_val->generators_G1[0]);

    element_t g1, g2;
    element_init_G1(g1, ret_val->pairing);
    element_init_G2(g2, ret_val->pairing);
    element_random(g1);
    element_random(g2);
    element_t gamma, omega;
    element_init_Zr(gamma, ret_val->pairing);
    element_init_G2(omega, ret_val->pairing);
    element_random(gamma);
    element_mul_zn(omega, g2, gamma);

    element_set(ret_val->G1, g1);
    element_set(ret_val->G2, g2);
    element_init_G2(ret_val->omega, ret_val->pairing);
    element_set(ret_val->omega, omega);
    element_init_Zr(sec_val->gamma, ret_val->pairing);
    element_set(sec_val->gamma, gamma);
    element_init_G1(ret_val->wf, ret_val->pairing);
    element_set(ret_val->wf, wf);
    element_printf("  Omega = %B\n", ret_val->omega);
    element_printf("  Water Function : %B\n", wf);
    long unsigned int a = mpz_get_ui(ret_val->a);
    mpz_t *S = (mpz_t *) malloc(sizeof(mpz_t) * a);
    element_t *gatt = (element_t *) malloc(sizeof(element_t) * a);
    printf("  Generator Attributes : For particlar attributes you can print\n");

    for (unsigned long int i = 0; i < a; i++) {
        mpz_init(S[i]);
        while (mpz_cmp_ui(S[i], 0) == 0)
            mpz_urandomm(S[i], state, ret_setup.p);
        element_init_G1(gatt[i], ret_setup.pairing);
        element_pow_mpz(gatt[i], g1, S[i]);
    }
    ret_val->gatt = gatt;
    sec_val->S = S;

    element_t alpha1, alpha2, t1, t2, alpha1_dash;
    element_init_Zr(alpha1, ret_val->pairing);
    element_init_Zr(alpha2, ret_val->pairing);
    element_init_Zr(alpha1_dash, ret_val->pairing);
    element_init_Zr(t1, ret_val->pairing);
    element_init_Zr(t2, ret_val->pairing);

    element_random(alpha2);
    element_random(alpha1);
    element_random(alpha1_dash);
    element_random(t1);
    element_random(t2);

    element_t u1_dash[2], v1[2], u1[2];
    element_t *u = (element_t *) malloc(sizeof(element_t) * 4);
    element_t *u_dash = (element_t *) malloc(sizeof(element_t) * 4);
    element_t *v = (element_t *) malloc(sizeof(element_t) * 4);

    for (int i = 0; i < 2; i++) {
        element_init_G1(u1[i], ret_setup.pairing);
        element_init_G1(u1_dash[i], ret_setup.pairing);
        element_init_G2(v1[i], ret_setup.pairing);
    }

    for (int i = 0; i < 4; i++) {
        element_init_G1(u[i], ret_setup.pairing);
        element_init_G1(u_dash[i], ret_setup.pairing);
        element_init_G2(v[i], ret_setup.pairing);
    }

    element_set(u1[0], g1);
    element_pow_zn(u1[1], g1, alpha1);
    element_set(u1_dash[0], g1);
    element_pow_zn(u1_dash[1], g1, alpha1_dash);
    element_set(v1[0], g2);
    element_pow_zn(v1[1], g2, alpha2);
    element_set(u[0], u1[0]);
    element_set(u[1], u1[1]);
    element_pow_zn(u[2], u1[0], t1);
    element_pow_zn(u[3], u1[1], t1);
    element_set(u_dash[0], u1_dash[0]);
    element_set(u_dash[1], u1_dash[1]);
    element_pow_zn(u_dash[2], u1_dash[0], t1);
    element_pow_zn(u_dash[3], u1_dash[1], t1);
    element_set(v[0], v1[0]);
    element_set(v[1], v1[1]);
    element_pow_zn(v[2], v1[0], t2);
    element_pow_zn(v[3], v1[1], t2);

    ret_setup.u = u;
    ret_setup.v = v;
    ret_setup.u_dash = u_dash;

    printf("  Groth Sahai Proof Vectors :\n");
    for (int i = 0; i < 4; i++) {
        element_printf("\tu[%d]: %B\n", i, u[i]);
        element_printf("\tu_dash[%d]: %B\n", i, u_dash[i]);
        element_printf("\tv[%d]: %B\n", i, v[i]);
    }

    element_printf("Issuing key : \n  Gamma = %B\n"
                   "  Set S = For particular attributes of secrets you can print\n",
                   gamma);
    element_printf("Opener key = %B\n", alpha1);
    element_printf("Tracer key = %B\n", alpha1_dash);
}

void join(public *retval, unsigned long int j, user_output *U_val, reg *rval, secret *sec_val) {
    printf("\n----------------Joint Algorithm for User U%lu------------------\n", j);

    U_val->i = j;
    rval->i = j;

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1));

    mpz_init(U_val->upki);
    mpz_init(rval->upki);
    mpz_init(U_val->uski);
    mpz_urandomm(U_val->upki, state, retval->p);
    mpz_set(rval->upki, U_val->upki);
    mpz_urandomm(U_val->uski, state, retval->p);

    element_t yi_dash, g1, Yi_dash;
    element_init_Zr(yi_dash, retval->pairing);
    element_init_G1(g1, retval->pairing);
    element_init_G1(Yi_dash, retval->pairing);
    element_random(yi_dash);
    element_set(g1, ret_setup.G1);
    element_mul_zn(Yi_dash, g1, yi_dash);
    element_printf("User U%lu sends :\n  Yi_dash : %B\n", j, Yi_dash);

    element_t xi, yi_dash_dash, Yi_dash_dash, Ai, gamma, add;
    element_init_Zr(xi, retval->pairing);
    element_init_Zr(yi_dash_dash, retval->pairing);
    element_init_Zr(U_val->yi_dash_dash, retval->pairing);
    element_init_G1(Yi_dash_dash, retval->pairing);
    element_init_G1(Ai, retval->pairing);
    element_init_Zr(gamma, retval->pairing);
    element_init_Zr(add, retval->pairing);
    element_random(xi);
    element_random(yi_dash_dash);
    element_set(U_val->yi_dash_dash, yi_dash_dash);
    element_pow_zn(Yi_dash_dash, g1, yi_dash_dash);
    element_set(gamma, sec_val->gamma);

    element_set(Ai, retval->generators_G1[0]);
    element_mul(Ai, Ai, Yi_dash);
    element_mul(Ai, Ai, Yi_dash_dash);
    element_add(add, gamma, xi);
    element_invert(add, add);
    element_pow_zn(Ai, Ai, add);

    element_t Xi_2;
    element_init_G2(Xi_2, retval->pairing);
    element_set(Xi_2, retval->G2);
    element_pow_zn(Xi_2, Xi_2, xi);

    mpz_t a;
    mpz_init(a);
    while (mpz_get_ui(a) == 0) {
        mpz_urandomm(a, state, retval->a);
    }

    long unsigned int *A = (long unsigned int *) malloc(sizeof(long unsigned int) * mpz_get_ui(a));
    long unsigned int k = 0;
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
    } while (k < mpz_get_ui(a));

    element_t temp_pow;
    element_init_Zr(temp_pow, retval->pairing);
    element_t *Ti = (element_t *) malloc(sizeof(element_t) * mpz_get_ui(retval->k));
    for (long unsigned int i = 0; i < mpz_get_ui(retval->a); i++) {
        element_init_G1(Ti[i], retval->pairing);
        element_mul_mpz(temp_pow, add, sec_val->S[i]);
        element_pow_zn(Ti[i], retval->generators_G1[0], temp_pow);
        element_set0(temp_pow);
    }

    element_printf("Group Manager sends :\n");
    element_printf("  y%lu_dash_dash : %B\n", j, yi_dash_dash);
    element_printf("  A%lu : %B\n", j, Ai);
    element_printf("  X%lu_2 : %B\n", j, Xi_2);
    for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
        long unsigned int x = A[i];
        element_printf("  T%lu_%lu : %B\n", j, x, Ti[x]);
    }

    element_t yi, Yi, temp1, temp2, temp3, temp4, temp_g2;
    element_init_Zr(yi, retval->pairing);
    element_init_G1(Yi, retval->pairing);
    element_init_GT(temp1, retval->pairing);
    element_init_GT(temp2, retval->pairing);
    element_init_GT(temp3, retval->pairing);
    element_init_GT(temp4, retval->pairing);
    element_init_G2(temp_g2, retval->pairing);
    element_add(yi, yi_dash, yi_dash_dash);
    element_mul(temp_g2, Xi_2, ret_setup.omega);
    element_pairing(temp1, Ai, temp_g2);

    element_pairing(temp2, retval->generators_G1[0], retval->G2);
    element_pairing(temp3, retval->G1, retval->G2);
    element_pow_zn(temp3, temp3, yi);
    element_mul(temp4, temp3, temp2);

    printf("User U%lu processing : \n", j);
    if (element_cmp(temp1, temp4) == 0)
        printf("\n   -----First Verification Done by user-----\n");
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");
    element_pow_zn(Yi, g1, yi);

    printf("\n  Digital Signature Sigma_%lu for : \n", j);
    element_printf("\tA%lu : %B\n", j, Ai);
    element_printf("\tX%lu_2 : %B\n", j, Xi_2);
    element_printf("\tY%lu : %B\n   Processing i.e. DSig_usk%lu :\n", j, Yi, j);

    mpz_t hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki;
    mpz_inits(hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki, NULL);
    element_to_mpz(hash_Ai, Ai);
    element_to_mpz(hash_Xi_2, Xi_2);
    element_to_mpz(hash_Yi, Yi);

    mpz_add(combined_hash, hash_Ai, hash_Xi_2);
    mpz_add(combined_hash, combined_hash, hash_Yi);
    mpz_add(combined_hash, combined_hash, U_val->upki);

    mpz_t n, d, e;
    mpz_inits(n, d, e, NULL);
    generate_rsa_keys(n, d, e, state);

    mpz_t signature;
    mpz_init(signature);
    mpz_init(U_val->signature);
    mpz_init(rval->signature);
    mpz_powm(signature, combined_hash, d, n);

    mpz_set(U_val->rsa1, d);
    mpz_set(U_val->rsa2, n);
    mpz_set(U_val->signature, signature);
    mpz_set(rval->signature, signature);

    mpz_t verified_message;
    mpz_init(verified_message);
    mpz_powm(verified_message, signature, e, n);
    mpz_set(rval->rsa2, n);
    mpz_set(rval->rsa3, e);

    printf("Group Manager GM processing : ");
    if (mpz_cmp(verified_message, combined_hash) == 0)
        printf("\n----GM verifies U%lu digital Signature----\n", j);
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");

    element_t Xi_1;
    element_init_G1(Xi_1, retval->pairing);
    element_pow_zn(Xi_1, g1, xi);
    element_printf("GM sends X%lu_1 = %B\n", j, Xi_1);

    printf("User U%lu processing : ", j);
    element_pairing(temp1, Xi_1, retval->G2);
    element_pairing(temp2, retval->G1, Xi_2);
    if (element_cmp(temp1, temp2) == 0)
        printf("\n----Verification passed : U%lu owns the membership----\n", j);
    else
        printf("\n-x-x-x-Verification failed-x-x-x-\n");

    printf("User%lu Secret Key : \n", j);
    printf("   Valid Membership Certificates : \n");
    element_printf("\tA%lu = %B\n", j, Ai);
    element_printf("\tX%lu = %B %B\n", j, Xi_1, Xi_2);
    element_printf("\ty%lu = %B\n", j, yi);
    printf("   Attribute certificate :\n");
    for (long unsigned int i = 0; i < mpz_get_ui(a); i++) {
        long unsigned int x = A[i];
        element_printf("  T%lu_%lu : %B\n", j, x, Ti[x]);
    }

    element_init_Zr(U_val->yi, retval->pairing);
    element_set(U_val->yi, yi);
    element_init_G1(U_val->Ai, retval->pairing);
    element_set(U_val->Ai, Ai);
    element_init_G2(U_val->Xi_1, retval->pairing);
    element_init_G2(U_val->Xi_2, retval->pairing);
    element_set(U_val->Xi_1, Xi_1);
    element_set(U_val->Xi_2, Xi_2);
    element_init_G1(U_val->Yi, retval->pairing);
    element_set(U_val->Yi, Yi);
    element_init_Zr(U_val->xi, retval->pairing);
    U_val->Ti = Ti;
    U_val->a = mpz_get_ui(a);
    U_val->A = A;
    element_set(U_val->xi, xi);

    element_init_G1(rval->Ai, retval->pairing);
    element_set(rval->Ai, Ai);
    element_init_G2(rval->Xi_1, retval->pairing);
    element_init_G2(rval->Xi_2, retval->pairing);
    element_set(rval->Xi_1, Xi_1);
    element_set(rval->Xi_2, Xi_2);
    element_init_G1(rval->Yi, retval->pairing);
    element_set(rval->Yi, Yi);

    printf("\n");
}

void buildtree(public *retval, mpz_t *S) {
    printf("\n----------------Build Tree Algorithm------------------\n");
    long unsigned int i = 0;
    mpz_t max, dummy, temp1;
    mpz_inits(max, dummy, temp1, NULL);
    while (1) {
        mpz_ui_pow_ui(temp1, 2, i);
        if (mpz_cmp(retval->a, temp1) > 0) {
            i++;
        } else
            break;
    }
    mpz_ui_pow_ui(dummy, 2, i);
    mpz_sub_ui(dummy, dummy, 1);
    mpz_ui_pow_ui(temp1, 2, i + 1);
    mpz_sub_ui(temp1, temp1, 1);
    mpz_add(max, temp1, dummy);

    Node *T = (Node *) malloc(mpz_get_ui(temp1) * sizeof(Node));
    Node *T_d = (Node *) malloc(mpz_get_ui(dummy) * sizeof(Node));
    for (i = 0; i < mpz_get_ui(temp1); i++) {
        initializeNode(&T[i]);
    }
    for (i = 0; i < mpz_get_ui(dummy); i++) {
        initializeNode(&T_d[i]);
    }

    for (long unsigned int j = mpz_get_ui(dummy); j < mpz_get_ui(temp1); j++) {
        i = j - mpz_get_ui(dummy);
        mpf_set_z(T[j].secret_value, S[i]);
    }

    mpz_set_ui(T[0].name, 1);
    long unsigned int k = 1, l = 0, m = 2;
    while (m <= mpz_get_ui(max)) {
        mpz_set_ui(T[k].name, m);
        mpz_set_ui(T[k + 1].name, m + 1);
        mpz_set_ui(T_d[l].name, m + 2);
        k = k + 2;
        l++;
        m = m + 3;
    }

    i = mpz_get_ui(temp1) - 1;
    while (i > 0) {
        long unsigned int j = (i / 2) - 1;
        mpf_t temp4, temp5;
        mpf_inits(temp4, temp5, NULL);
        mpf_set_z(temp4, T[i - 1].name);
        mpf_set_z(temp5, T[i].name);
        mpf_mul(temp4, T[i].secret_value, temp4);
        mpf_mul(temp5, T[i - 1].secret_value, temp5);
        mpf_sub(T[j].secret_value, temp5, temp4);

        mpf_mul_ui(temp4, T[i].secret_value, 2);
        mpf_sub(T_d[j].secret_value, temp4, T[i - 1].secret_value);
        i = i - 2;
    }
    l = 0, m = 1, k = 1;
    gmp_printf("Root  : %.Ff\n", T[0].secret_value);
    while (m < mpz_get_ui(max)) {
        gmp_printf("Node %Zd: %.Ff\n", T[k].name, T[k].secret_value);
        gmp_printf("Node %Zd: %.Ff\n", T[k + 1].name, T[k + 1].secret_value);
        k = k + 2;
        gmp_printf("Node %Zd: %.Ff\n", T_d[l].name, T_d[l].secret_value);
        l++;
        m = m + 3;
    }
    mpz_init_set(retval->max, max);
    mpz_init_set(retval->temp1, temp1);
    mpz_init_set(retval->d, dummy);
    retval->T_d = T_d;
    private.T = T;
    element_init_G2(retval->vT, retval->pairing);
    mpz_set_f(temp1, T[0].secret_value);
    mpz_abs(temp1, temp1);
    element_t tem1;
    element_init_Zr(tem1, ret_setup.pairing);
    element_set_mpz(tem1, temp1);
    element_pow_zn(retval->vT, retval->G2, tem1);
    printf("\n");
}

void buildtreevalidity(public *retval) {
    printf("\n----------------Build Tree Validity Algorithm------------------\n");
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (rand() + 1) * (rand() + 1));
    mpz_t a;
    mpz_init(a);
    while (mpz_get_ui(a) == 0) {
        mpz_urandomm(a, state, retval->a);
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
    } while (k < mpz_get_ui(a));
    insertionSort(A, mpz_get_ui(a));
    for (k = 0; k < mpz_get_ui(a); k++) {
        gmp_printf("\t%lu", A[k] + 1);
    }
    verify_tree(A, retval->d, retval->temp1, retval->max, ret_setup.T_d, 0);
}

void sign(user_output *U_val) {
    gmp_printf("\n*----------------Sign Algorithm for U%lu------------------\n", U_val->i);
    printf("Attributes :");
    for (long unsigned int i = 0; i < U_val->a; i++) {
        gmp_printf("\t%lu", U_val->A[i] + 1);
    }
    verify_tree(U_val->A, ret_setup.d, ret_setup.temp1, ret_setup.max, ret_setup.T_d, U_val->i);

    element_t temp, ID, z, r;
    element_init_G1(ID, ret_setup.pairing);
    element_init_Zr(temp, ret_setup.pairing);
    element_init_Zr(z, ret_setup.pairing);
    element_init_Zr(r, ret_setup.pairing);
    element_init_Zr(R[U_val->i - 1].z, ret_setup.pairing);
    element_init_Zr(R[U_val->i - 1].r, ret_setup.pairing);

    element_random(z);
    element_add(temp, z, U_val->yi);
    element_invert(temp, temp);
    element_pow_zn(ID, ret_setup.G1, temp);
    element_random(r);
    element_set(R[U_val->i - 1].r, r);
    element_set(R[U_val->i - 1].z, z);

    element_t *rho3 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *rho7 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t rho4, rho5, rho6, rho8, rho9, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, vT;
    element_init_G1(rho3[0], ret_setup.pairing);
    element_init_G2(rho3[1], ret_setup.pairing);
    element_init_G1(rho4, ret_setup.pairing);
    element_init_G1(rho5, ret_setup.pairing);
    element_init_G1(rho6, ret_setup.pairing);
    element_init_G1(rho7[0], ret_setup.pairing);
    element_init_G2(rho7[1], ret_setup.pairing);
    element_init_G1(rho8, ret_setup.pairing);
    element_init_G2(rho9, ret_setup.pairing);
    element_init_Zr(temp1, ret_setup.pairing);
    element_init_Zr(temp2, ret_setup.pairing);
    element_init_G1(temp3, ret_setup.pairing);
    element_init_G1(temp4, ret_setup.pairing);
    element_init_G2(temp5, ret_setup.pairing);
    element_init_G2(temp9, ret_setup.pairing);
    element_init_GT(temp6, ret_setup.pairing);
    element_init_GT(temp7, ret_setup.pairing);
    element_init_GT(temp8, ret_setup.pairing);
    element_init_G2(vT, ret_setup.pairing);
    element_init_G2(U_val->vT, ret_setup.pairing);

    element_pow_mpz(vT, ret_setup.G2, R[U_val->i - 1].sT);
    element_set(U_val->vT, vT);

    element_set(rho3[0], U_val->Xi_1);
    element_set(rho3[1], U_val->Xi_2);

    element_add(temp1, private.gamma, U_val->xi);
    element_invert(temp2, temp1);
    if (mpz_cmp_ui(R[U_val->i - 1].sT1, 0) > 0) {
        element_mul_mpz(temp1, temp2, R[U_val->i - 1].sT1);
        element_pow_zn(rho4, ret_setup.generators_G1[0], temp1);
    } else {
        mpz_t tem1;
        mpz_init(tem1);
        mpz_abs(tem1, R[U_val->i - 1].sT1);
        element_mul_mpz(temp2, temp2, tem1);
        element_pow_zn(rho4, ret_setup.generators_G1[0], temp2);
        element_invert(rho4, rho4);
    }

    if (mpz_cmp_ui(R[U_val->i - 1].sT2, 0) > 0) {
        element_pow_mpz(rho5, ret_setup.generators_G1[0], R[U_val->i - 1].sT2);
    } else {
        mpz_t tem1;
        mpz_init_set(tem1, R[U_val->i - 1].sT2);
        mpz_abs(tem1, tem1);
        element_pow_mpz(rho5, ret_setup.generators_G1[0], tem1);
        element_invert(rho5, rho5);
    }

    element_set(rho6, ID);

    element_pow_zn(rho7[0], ret_setup.G1, z);
    element_pow_zn(rho7[1], ret_setup.G2, z);

    element_set(temp3, ret_setup.wf);
    element_pow_zn(temp3, temp3, r);
    element_pow_zn(temp4, ret_setup.generators_G1[0], z);
    element_mul(rho8, temp3, temp4);

    element_pow_zn(rho9, ret_setup.G2, r);

    element_mul(temp5, ret_setup.omega, rho3[1]);
    element_pairing(temp6, U_val->Ai, temp5);
    element_pairing(temp7, ret_setup.generators_G1[0], ret_setup.G2);
    element_pow_zn(temp9, ret_setup.G2, U_val->yi);
    element_pairing(temp8, ret_setup.G1, temp9);
    element_mul(temp7, temp7, temp8);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----Signer has valid Membership Cerficate-----");
    else
        printf("-*-*-Fraud-*-*-");

    element_pairing(temp6, rho4, temp5);
    element_pairing(temp7, rho5, ret_setup.G2);
    element_pairing(temp8, ret_setup.generators_G1[0], U_val->vT);
    element_mul(temp6, temp6, temp7);
    if (element_cmp(temp6, temp8) == 0)
        printf("\n-----Signer has valid Attribute Cerficate-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pow_zn(temp5, ret_setup.G2, U_val->yi);
    element_set1(temp9);
    element_mul(temp9, temp5, rho7[1]);
    element_pairing(temp6, rho6, temp9);
    element_pairing(temp7, ret_setup.G1, ret_setup.G2);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----Signer has valid ID-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho8, ret_setup.G2);
    element_pairing(temp7, ret_setup.generators_G1[0], rho7[1]);
    element_pairing(temp8, ret_setup.wf, rho9);
    element_mul(temp7, temp7, temp8);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----Water Signature valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pow_zn(temp4, ret_setup.G1, U_val->yi);
    element_pairing(temp6, temp4, ret_setup.G2);
    element_pow_zn(temp5, ret_setup.G2, U_val->yi);
    element_pairing(temp7, ret_setup.G1, temp5);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----yi commited valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho3[0], ret_setup.G2);
    element_pairing(temp7, ret_setup.G1, rho3[1]);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----Xi commited valid-----");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_pairing(temp6, rho7[0], ret_setup.G2);
    element_pairing(temp7, ret_setup.G1, rho7[1]);
    if (element_cmp(temp6, temp7) == 0)
        printf("\n-----Non-frameability Adversary checked valid-----\n");
    else
        printf("\n-*-*-Fraud-*-*-");

    element_init_G1(U_val->rho1, ret_setup.pairing);
    element_init_Zr(U_val->rho2, ret_setup.pairing);
    element_init_G1(U_val->rho4, ret_setup.pairing);
    element_init_G1(U_val->rho5, ret_setup.pairing);
    element_init_G1(U_val->rho6, ret_setup.pairing);
    element_init_G1(U_val->rho8, ret_setup.pairing);
    element_init_G2(U_val->rho9, ret_setup.pairing);
    element_set(U_val->rho1, U_val->Ai);
    element_set(U_val->rho2, U_val->yi);
    element_set(U_val->rho4, rho4);
    element_set(U_val->rho5, rho5);
    element_set(U_val->rho6, rho6);
    element_set(U_val->rho8, rho8);
    element_set(U_val->rho9, rho9);
    U_val->rho3 = rho3;
    U_val->rho7 = rho7;

    element_printf("\n  User U%lu Signature :", U_val->i);
    element_printf("\n\trho1 : %B", U_val->Ai);
    element_printf("\n\trho2 : %B", U_val->yi);
    element_printf("\n\trho3,1 : %B", rho3[0]);
    element_printf("\n\trho3,2 : %B", rho3[1]);
    element_printf("\n\trho4 : %B", rho4);
    element_printf("\n\trho5 : %B", rho5);
    element_printf("\n\trho6 : %B", rho6);
    element_printf("\n\trho7,1 : %B", rho7[0]);
    element_printf("\n\trho7,2 : %B", rho7[1]);
    element_printf("\n\trho8 : %B", rho8);
    element_printf("\n\trho9 : %B", rho9);

    element_t *sigma1 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma2_1 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma2_2 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma3_1 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma3_2 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma4 = (element_t *) malloc(sizeof(element_t) * 2);
    element_t *sigma5 = (element_t *) malloc(sizeof(element_t) * 2);

    element_t tem1[2], tem2[2], tem3[2], tem4[2], tem5[2], tem6[2], tem7[2], tem8[2];
    for (int i = 0; i < 2; i++) {
        element_init_G1(tem1[i], ret_setup.pairing);
        element_init_G1(tem2[i], ret_setup.pairing);
        element_init_G1(tem3[i], ret_setup.pairing);
        element_init_G1(tem4[i], ret_setup.pairing);
        element_init_G2(tem5[i], ret_setup.pairing);
        element_init_G2(tem6[i], ret_setup.pairing);
        element_init_G2(tem7[i], ret_setup.pairing);
        element_init_G1(sigma1[i], ret_setup.pairing);
        element_init_G1(sigma2_1[i], ret_setup.pairing);
        element_init_G2(sigma2_2[i], ret_setup.pairing);
        element_init_G1(sigma3_1[i], ret_setup.pairing);
        element_init_G2(sigma3_2[i], ret_setup.pairing);
        element_init_G1(sigma4[i], ret_setup.pairing);
        element_init_G1(sigma5[i], ret_setup.pairing);
    }

    element_init_G1(tem8[0], ret_setup.pairing);
    element_init_G2(tem8[1], ret_setup.pairing);
    element_set_si(tem8[0], 1);
    element_set_si(tem8[1], 1);

    element_pow_zn(tem1[0], ret_setup.u[0], r);
    element_pow_zn(tem1[1], ret_setup.u[1], r);
    element_pow_zn(tem2[0], ret_setup.u[2], z);
    element_pow_zn(tem2[1], ret_setup.u[3], z);
    element_mul(tem4[0], tem2[0], tem1[0]);
    element_mul(tem4[1], tem2[1], tem1[1]);
    element_mul(sigma1[0], tem4[0], tem8[0]);
    element_mul(sigma1[1], tem4[1], U_val->rho1);
    element_printf("\n\tsigma1 : %B %B", sigma1[0], sigma1[1]);

    element_mul(tem2[0], ret_setup.u[2], tem8[0]);
    element_mul(tem2[1], ret_setup.u[3], ret_setup.G1);
    element_pow_zn(tem3[0], tem2[0], U_val->rho2);
    element_pow_zn(tem3[1], tem2[1], U_val->rho2);
    element_mul(sigma2_1[0], tem3[0], tem1[0]);
    element_mul(sigma2_1[1], tem3[1], tem2[1]);
    element_printf("\n\tsigma2_1 : %B %B", sigma2_1[0], sigma2_1[1]);

    element_pow_zn(tem5[0], ret_setup.v[0], r);
    element_pow_zn(tem5[1], ret_setup.v[1], r);
    element_mul(tem6[0], ret_setup.v[2], tem8[1]);
    element_mul(tem6[1], ret_setup.v[3], ret_setup.G2);
    element_pow_zn(tem7[0], tem6[0], U_val->rho2);
    element_pow_zn(tem7[1], tem6[1], U_val->rho2);
    element_mul(sigma2_2[0], tem7[0], tem5[0]);
    element_mul(sigma2_2[1], tem7[1], tem5[1]);
    element_printf("\n\tsigma2_2 : %B %B", sigma2_2[0], sigma2_2[1]);

    element_mul(sigma3_1[0], tem4[0], tem8[0]);
    element_mul(sigma3_1[1], tem4[1], U_val->rho3[0]);
    element_printf("\n\tsigma3,1 : %B %B", sigma3_1[0], sigma3_1[1]);
    element_pow_zn(tem5[0], ret_setup.v[0], r);
    element_pow_zn(tem5[1], ret_setup.v[1], r);
    element_pow_zn(tem6[0], ret_setup.v[2], z);
    element_pow_zn(tem6[1], ret_setup.v[3], z);
    element_mul(tem6[0], tem6[0], tem5[0]);
    element_mul(tem6[1], tem6[1], tem5[1]);
    element_mul(sigma3_2[0], tem6[0], tem8[1]);
    element_mul(sigma3_2[1], tem6[1], U_val->rho3[1]);
    element_printf("\n\tsigma3,2 : %B %B", sigma3_2[0], sigma3_2[1]);

    element_pow_zn(tem1[0], ret_setup.u_dash[0], r);
    element_pow_zn(tem1[1], ret_setup.u_dash[1], r);
    element_pow_zn(tem2[0], ret_setup.u_dash[2], z);
    element_pow_zn(tem2[1], ret_setup.u_dash[3], z);
    element_mul(tem4[0], tem2[0], tem1[0]);
    element_mul(tem4[1], tem2[1], tem1[1]);
    element_mul(sigma4[0], tem4[0], tem8[0]);
    element_mul(sigma4[1], tem4[1], U_val->rho4);
    element_printf("\n\tsigma4 : %B %B", sigma4[0], sigma4[1]);

    element_mul(sigma5[0], tem4[0], tem8[0]);
    element_mul(sigma5[1], tem4[1], U_val->rho5);
    element_printf("\n\tsigma5 : %B %B\n", sigma5[0], sigma5[1]);

    U_val->sigma1 = sigma1;
    U_val->sigma2_1 = sigma2_1;
    U_val->sigma2_2 = sigma2_2;
    U_val->sigma3_1 = sigma3_1;
    U_val->sigma3_2 = sigma3_2;
    U_val->sigma4 = sigma4;
    U_val->sigma5 = sigma5;
}

void verify(user_output *U_val) {
    mpz_t temp1;
    mpz_init(temp1);
    element_t groot, vT, tem1, tem2, sT;
    element_init_G1(groot, ret_setup.pairing);
    element_init_G2(vT, ret_setup.pairing);
    element_init_GT(tem1, ret_setup.pairing);
    element_init_GT(tem2, ret_setup.pairing);
    element_init_Zr(sT, ret_setup.pairing);
    mpz_set_f(temp1, U_val->T[0].secret_value);
    mpz_abs(temp1, temp1);
    element_set_mpz(sT, temp1);
    element_pow_zn(groot, ret_setup.G1, sT);
    element_set(vT, ret_setup.vT);
    element_pairing(tem1, groot, ret_setup.G2);
    element_pairing(tem2, ret_setup.G1, vT);
    if (element_cmp(tem1, tem2) == 0)
        gmp_printf("\n-----Verify algorithm valid for User%lu-----\n", U_val->i);

    else
        printf("\nx-x-x-Verify algorithm invalid-x-x-x\n");
}

void openuser(user_output *U_val, reg *rval, element_t ok_user) {
    if (private.alpha1 == ok_user) {
        printf("\n*----------OpenUser verified-----------");

        mpz_t hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki;
        mpz_inits(hash_Ai, hash_Xi_2, hash_Yi, combined_hash, hash_upki, NULL);
        element_to_mpz(hash_Ai, U_val->Ai);
        element_to_mpz(hash_Xi_2, U_val->Xi_2);
        element_to_mpz(hash_Yi, U_val->Yi);

        mpz_add(combined_hash, hash_Ai, hash_Xi_2);
        mpz_add(combined_hash, combined_hash, hash_Yi);
        mpz_add(combined_hash, combined_hash, rval->upki);
        mpz_t verified_message;
        mpz_init(verified_message);
        mpz_powm(verified_message, rval->signature, rval->rsa3, rval->rsa2);
        if (mpz_cmp(verified_message, combined_hash) == 0) {
            printf("\n  ----OpenUser verifies U%lu digital Signature----\n", rval->i);

            element_t *sigma1_verify = (element_t *) malloc(sizeof(element_t) * 2);
            element_t tem1[2], tem2[2], tem4[2];
            for (int i = 0; i < 2; i++) {
                element_init_G1(tem1[i], ret_setup.pairing);
                element_init_G1(tem2[i], ret_setup.pairing);
                element_init_G1(tem4[i], ret_setup.pairing);
                element_init_G1(sigma1_verify[i], ret_setup.pairing);
            }
            element_pow_zn(tem1[0], ret_setup.u[0], rval->r);
            element_pow_zn(tem1[1], ret_setup.u[1], rval->r);
            element_pow_zn(tem2[0], ret_setup.u[2], rval->z);
            element_pow_zn(tem2[1], ret_setup.u[3], rval->z);
            element_mul(tem4[0], tem2[0], tem1[0]);
            element_mul(tem4[1], tem2[1], tem1[1]);
            element_div(sigma1_verify[0], U_val->sigma1[0], tem4[0]);
            element_div(sigma1_verify[1], U_val->sigma1[1], tem4[1]);

            if (element_cmp(rval->Ai, sigma1_verify[1]) == 0) {
                gmp_printf("\t---User Identity %lu ---\n", rval->i);
            } else
                printf("\n\t-*-*-Wrong or No User Identity-*-*- ");
        } else
            printf("\n-x-x-x-Verification failed-x-x-x-\n");
    } else
        printf("\n-x-x-x-Open User Verification failed-x-x-x-\n");
}

void traceatt(user_output *U_val, reg *rval, element_t tk_att) {
    if (private.alpha1_dash == tk_att) {
        gmp_printf("\n----Tracer is valid----\n");

        element_t *sigma5_verify = (element_t *) malloc(sizeof(element_t) * 2);
        element_t tem1[2], tem2[2], tem4[2];
        for (int i = 0; i < 2; i++) {
            element_init_G1(tem1[i], ret_setup.pairing);
            element_init_G1(tem2[i], ret_setup.pairing);
            element_init_G1(tem4[i], ret_setup.pairing);
            element_init_G1(sigma5_verify[i], ret_setup.pairing);
        }
        element_pow_zn(tem1[0], ret_setup.u_dash[0], rval->r);
        element_pow_zn(tem1[1], ret_setup.u_dash[1], rval->r);
        element_pow_zn(tem2[0], ret_setup.u_dash[2], rval->z);
        element_pow_zn(tem2[1], ret_setup.u_dash[3], rval->z);
        element_mul(tem4[0], tem2[0], tem1[0]);
        element_mul(tem4[1], tem2[1], tem1[1]);
        element_div(sigma5_verify[0], U_val->sigma5[0], tem4[0]);
        element_div(sigma5_verify[1], U_val->sigma5[1], tem4[1]);

        element_t rho5_verify;
        element_init_G1(rho5_verify, ret_setup.pairing);
        if (mpz_cmp_ui(rval[U_val->i - 1].sT2, 0) > 0) {
            element_pow_mpz(rho5_verify, ret_setup.generators_G1[0], R[U_val->i - 1].sT2);
        } else {
            mpz_t temp1;
            mpz_init_set(temp1, R[U_val->i - 1].sT2);
            mpz_abs(temp1, temp1);
            element_pow_mpz(rho5_verify, ret_setup.generators_G1[0], temp1);
            element_invert(rho5_verify, rho5_verify);
        }

        gmp_printf("U%lu Attribute Name: ", U_val->i);
        if (element_cmp(rho5_verify, sigma5_verify[1]) == 0) {
            for (long unsigned int i = 0; i < (U_val->a); i++) {
                gmp_printf("   %lu", U_val->A[i] + 1);
            }
            printf("\n");
        } else
            printf("\n\t-*-*-Wrong or No User Identity-*-*- ");
    } else
        gmp_printf("\n----Tracer is invalid----\n");
}

int main() {
    srand(time(NULL));
    mpz_t security_parameter;
    mpz_init(security_parameter);
    mpz_set_ui(security_parameter, 16);

    setup(&ret_setup, security_parameter);

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
