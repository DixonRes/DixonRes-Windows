#include "fq_nmod_roots.h"

// nmod_roots structure management

void nmod_roots_init(nmod_roots_t roots) {
    roots->alloc = 4;
    roots->roots = (mp_limb_t*)malloc(roots->alloc * sizeof(mp_limb_t));
    roots->mult = (slong*)malloc(roots->alloc * sizeof(slong));
    roots->num = 0;
}

void nmod_roots_clear(nmod_roots_t roots) {
    free(roots->roots);
    free(roots->mult);
    roots->num = 0;
    roots->alloc = 0;
}

void nmod_roots_fit_length(nmod_roots_t roots, slong len) {
    if (len > roots->alloc) {
        roots->alloc = FLINT_MAX(len, 2 * roots->alloc);
        roots->roots = (mp_limb_t*)realloc(roots->roots, roots->alloc * sizeof(mp_limb_t));
        roots->mult = (slong*)realloc(roots->mult, roots->alloc * sizeof(slong));
    }
}

void nmod_roots_add(nmod_roots_t roots, mp_limb_t root, slong mult) {
    nmod_roots_fit_length(roots, roots->num + 1);
    roots->roots[roots->num] = root;
    roots->mult[roots->num] = mult;
    roots->num++;
}

// fq_nmod_roots structure management

void fq_nmod_roots_init(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx) {
    roots->alloc = 4;
    roots->roots = (fq_nmod_struct*)malloc(roots->alloc * sizeof(fq_nmod_struct));
    roots->mult = (slong*)malloc(roots->alloc * sizeof(slong));
    roots->num = 0;
    roots->p = fq_nmod_ctx_prime(ctx);
    roots->d = fq_nmod_ctx_degree(ctx);
    
    // Initialize root elements
    for (int i = 0; i < roots->alloc; i++) {
        fq_nmod_init(roots->roots + i, ctx);
    }
}

void fq_nmod_roots_clear(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx) {
    for (int i = 0; i < roots->alloc; i++) {
        fq_nmod_clear(roots->roots + i, ctx);
    }
    free(roots->roots);
    free(roots->mult);
    roots->num = 0;
    roots->alloc = 0;
}

void fq_nmod_roots_fit_length(fq_nmod_roots_t roots, slong len, const fq_nmod_ctx_t ctx) {
    if (len > roots->alloc) {
        slong old_alloc = roots->alloc;
        roots->alloc = FLINT_MAX(len, 2 * roots->alloc);
        roots->roots = (fq_nmod_struct*)realloc(roots->roots, roots->alloc * sizeof(fq_nmod_struct));
        roots->mult = (slong*)realloc(roots->mult, roots->alloc * sizeof(slong));
        
        // Initialize new elements
        for (slong i = old_alloc; i < roots->alloc; i++) {
            fq_nmod_init(roots->roots + i, ctx);
        }
    }
}

void fq_nmod_roots_add(fq_nmod_roots_t roots, const fq_nmod_t root, slong mult, const fq_nmod_ctx_t ctx) {
    fq_nmod_roots_fit_length(roots, roots->num + 1, ctx);
    fq_nmod_set(roots->roots + roots->num, root, ctx);
    roots->mult[roots->num] = mult;
    roots->num++;
}

// Utility functions

double get_time_roots() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// nmod_poly version implementations

void nmod_simple_power_x_mod(nmod_poly_t result, mp_limb_t exp, const nmod_poly_t modulus) {
    if (exp == 0) {
        nmod_poly_one(result);
        return;
    }
    
    nmod_poly_t base;
    nmod_poly_init(base, modulus->mod.n);
    
    nmod_poly_set_coeff_ui(base, 1, 1);  // base = x
    nmod_poly_one(result);  // result = 1
    
    mp_limb_t e = exp;
    while (e > 0) {
        if (e & 1) {
            nmod_poly_mulmod(result, result, base, modulus);
        }
        e >>= 1;
        if (e > 0) {
            nmod_poly_mulmod(base, base, base, modulus);
        }
    }
    
    nmod_poly_clear(base);
}

void nmod_extract_linear_factors(nmod_roots_t roots, const nmod_poly_t poly, flint_rand_t state) {
    slong deg = nmod_poly_degree(poly);
    
    if (deg <= 0) return;
    
    if (deg == 1) {
        mp_limb_t a = nmod_poly_get_coeff_ui(poly, 1);
        mp_limb_t b = nmod_poly_get_coeff_ui(poly, 0);
        
        if (a != 0) {
            mp_limb_t inv_a = n_invmod(a, poly->mod.n);
            mp_limb_t root = n_mulmod2_preinv(poly->mod.n - b, inv_a, poly->mod.n, poly->mod.ninv);
            nmod_roots_add(roots, root, 1);
        }
        return;
    }
    
    // Randomized splitting
    nmod_poly_t g, h, r;
    nmod_poly_init(g, poly->mod.n);
    nmod_poly_init(h, poly->mod.n);
    nmod_poly_init(r, poly->mod.n);
    
    for (int attempt = 0; attempt < 10; attempt++) {
        nmod_poly_randtest(g, state, deg);
        
        mp_limb_t exp = (poly->mod.n - 1) / 2;
        nmod_poly_powmod_ui_binexp(h, g, exp, poly);
        
        mp_limb_t constant = nmod_poly_get_coeff_ui(h, 0);
        nmod_poly_set_coeff_ui(h, 0, n_submod(constant, 1, poly->mod.n));
        
        nmod_poly_gcd(r, h, poly);
        slong r_deg = nmod_poly_degree(r);
        
        if (r_deg > 0 && r_deg < deg) {
            nmod_extract_linear_factors(roots, r, state);
            
            nmod_poly_t quotient;
            nmod_poly_init(quotient, poly->mod.n);
            nmod_poly_div(quotient, poly, r);
            nmod_extract_linear_factors(roots, quotient, state);
            nmod_poly_clear(quotient);
            break;
        }
    }
    
    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(r);
}

slong our_nmod_poly_roots(nmod_roots_t roots, const nmod_poly_t poly, int with_multiplicity) {
    slong deg = nmod_poly_degree(poly);
    
    if (deg <= 0) return 0;
    
    if (deg == 1) {
        mp_limb_t a = nmod_poly_get_coeff_ui(poly, 1);
        mp_limb_t b = nmod_poly_get_coeff_ui(poly, 0);
        
        if (a != 0) {
            mp_limb_t inv_a = n_invmod(a, poly->mod.n);
            mp_limb_t root = n_mulmod2_preinv(poly->mod.n - b, inv_a, poly->mod.n, poly->mod.ninv);
            nmod_roots_add(roots, root, 1);
            return 1;
        }
        return 0;
    }
    
    // Cantor-Zassenhaus algorithm core
    nmod_poly_t x_to_p, x, frobenius, root_poly;
    nmod_poly_init(x_to_p, poly->mod.n);
    nmod_poly_init(x, poly->mod.n);
    nmod_poly_init(frobenius, poly->mod.n);
    nmod_poly_init(root_poly, poly->mod.n);
    
    nmod_simple_power_x_mod(x_to_p, poly->mod.n, poly);
    nmod_poly_set_coeff_ui(x, 1, 1);
    nmod_poly_sub(frobenius, x_to_p, x);
    nmod_poly_gcd(root_poly, poly, frobenius);
    
    if (nmod_poly_degree(root_poly) > 0) {
        flint_rand_t state;
        flint_rand_init(state);
        nmod_extract_linear_factors(roots, root_poly, state);
        flint_rand_clear(state);
    }
    
    nmod_poly_clear(x_to_p);
    nmod_poly_clear(x);
    nmod_poly_clear(frobenius);
    nmod_poly_clear(root_poly);
    
    return roots->num;
}

void fq_nmod_simple_power_x_mod(fq_nmod_poly_t result, mp_limb_t exp, const fq_nmod_poly_t modulus, const fq_nmod_ctx_t ctx) {
    if (exp == 0) {
        fq_nmod_poly_one(result, ctx);
        return;
    }
    
    fq_nmod_poly_t base;
    fq_nmod_poly_init(base, ctx);
    fq_nmod_poly_gen(base, ctx);  // base = x
    
    fq_nmod_poly_one(result, ctx);
    
    while (exp > 0) {
        if (exp & 1) {
            fq_nmod_poly_mulmod(result, result, base, modulus, ctx);
        }
        exp >>= 1;
        if (exp > 0) {
            fq_nmod_poly_mulmod(base, base, base, modulus, ctx);
        }
    }
    
    fq_nmod_poly_clear(base, ctx);
}

// Also fix the extract_linear_factors function
void fq_nmod_extract_linear_factors(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, flint_rand_t state, const fq_nmod_ctx_t ctx) {
    slong deg = fq_nmod_poly_degree(poly, ctx);
    
    if (deg <= 0) return;
    
    if (deg == 1) {
        fq_nmod_t a, b, root;
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fq_nmod_init(root, ctx);
        
        fq_nmod_poly_get_coeff(a, poly, 1, ctx);
        fq_nmod_poly_get_coeff(b, poly, 0, ctx);
        
        if (!fq_nmod_is_zero(a, ctx)) {
            fq_nmod_inv(root, a, ctx);
            fq_nmod_neg(b, b, ctx);
            fq_nmod_mul(root, root, b, ctx);
            fq_nmod_roots_add(roots, root, 1, ctx);
        }
        
        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(root, ctx);
        return;
    }
    
    // Random splitting with correct field size
    fq_nmod_poly_t g, h, r;
    fq_nmod_poly_init(g, ctx);
    fq_nmod_poly_init(h, ctx);
    fq_nmod_poly_init(r, ctx);
    
    // Get field order q = p^d
    fmpz_t q, exp;
    fmpz_init(q);
    fmpz_init(exp);
    fq_nmod_ctx_order(q, ctx);
    fmpz_sub_ui(exp, q, 1);
    fmpz_fdiv_q_2exp(exp, exp, 1);  // (q-1)/2
    
    int success = 0;
    for (int attempt = 0; attempt < 20 && !success; attempt++) {
        fq_nmod_poly_randtest(g, state, deg, ctx);
        
        fq_nmod_poly_t one_poly;
        fq_nmod_poly_init(one_poly, ctx);
        fq_nmod_poly_one(one_poly, ctx);
        
        // Compute g^((q-1)/2) - 1
        fq_nmod_poly_powmod_fmpz_binexp(h, g, exp, poly, ctx);
        fq_nmod_poly_sub(h, h, one_poly, ctx);
        
        fq_nmod_poly_gcd(r, h, poly, ctx);
        slong r_deg = fq_nmod_poly_degree(r, ctx);
        
        if (r_deg > 0 && r_deg < deg) {
            fq_nmod_extract_linear_factors(roots, r, state, ctx);
            
            fq_nmod_poly_t quotient;
            fq_nmod_poly_init(quotient, ctx);
            fq_nmod_poly_div(quotient, poly, r, ctx);
            fq_nmod_extract_linear_factors(roots, quotient, state, ctx);
            fq_nmod_poly_clear(quotient, ctx);
            success = 1;
        }
        
        fq_nmod_poly_clear(one_poly, ctx);
    }
    
    fmpz_clear(q);
    fmpz_clear(exp);
    fq_nmod_poly_clear(g, ctx);
    fq_nmod_poly_clear(h, ctx);
    fq_nmod_poly_clear(r, ctx);
}

slong our_fq_nmod_poly_roots(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, int with_multiplicity, const fq_nmod_ctx_t ctx) {
    slong deg = fq_nmod_poly_degree(poly, ctx);
    //printf("Input polynomial degree: %ld\n", deg);
    
    if (deg <= 0) return 0;
    
    if (deg < 1000) {
        fq_nmod_poly_factor_t factors;
        fq_nmod_poly_factor_init(factors, ctx);
        
        fq_nmod_poly_roots(factors, poly, with_multiplicity, ctx);
        
        for (slong i = 0; i < factors->num; i++) {
            slong factor_deg = fq_nmod_poly_degree(factors->poly + i, ctx);
            
            if (factor_deg == 1) {
                fq_nmod_t a, b, root;
                fq_nmod_init(a, ctx);
                fq_nmod_init(b, ctx);
                fq_nmod_init(root, ctx);
                
                fq_nmod_poly_get_coeff(a, factors->poly + i, 1, ctx);
                fq_nmod_poly_get_coeff(b, factors->poly + i, 0, ctx);
                
                if (!fq_nmod_is_zero(a, ctx)) {
                    fq_nmod_inv(root, a, ctx);
                    fq_nmod_neg(b, b, ctx);
                    fq_nmod_mul(root, root, b, ctx);
                    fq_nmod_roots_add(roots, root, with_multiplicity ? factors->exp[i] : 1, ctx);
                }
                
                fq_nmod_clear(a, ctx);
                fq_nmod_clear(b, ctx);
                fq_nmod_clear(root, ctx);
            }
        }
        
        fq_nmod_poly_factor_clear(factors, ctx);
        return roots->num;
    }
    
    if (deg == 1) {
        //printf("Linear case detected\n");
        fq_nmod_t a, b, root;
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx); 
        fq_nmod_init(root, ctx);
        
        fq_nmod_poly_get_coeff(a, poly, 1, ctx);
        fq_nmod_poly_get_coeff(b, poly, 0, ctx);
        
        if (!fq_nmod_is_zero(a, ctx)) {
            fq_nmod_inv(root, a, ctx);
            fq_nmod_neg(b, b, ctx);
            fq_nmod_mul(root, root, b, ctx);
            fq_nmod_roots_add(roots, root, 1, ctx);
            //printf("Found linear root: ");
            //fq_nmod_print_pretty(root, ctx);
            //printf("\n");
        }
        
        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(root, ctx);
        return roots->num;
    }
    
    //printf("Starting Cantor-Zassenhaus algorithm\n");
    
    mp_limb_t p = fq_nmod_ctx_prime(ctx);
    slong d = fq_nmod_ctx_degree(ctx);
    //printf("Field: F_{%lu^%ld}\n", p, d);
    
    fq_nmod_poly_t x_to_q, x, frobenius, splitting_poly;
    fq_nmod_poly_init(x_to_q, ctx);
    fq_nmod_poly_init(x, ctx);
    fq_nmod_poly_init(frobenius, ctx);
    fq_nmod_poly_init(splitting_poly, ctx);
    
    fmpz_t q;
    fmpz_init(q);
    fq_nmod_ctx_order(q, ctx);
    //printf("Field size q = ");
    //fmpz_print(q);
    //printf("\n");
    
    fq_nmod_poly_gen(x, ctx);  // x
    fq_nmod_poly_powmod_fmpz_binexp(x_to_q, x, q, poly, ctx);
    
    fq_nmod_poly_sub(frobenius, x_to_q, x, ctx);
    
    //printf("Frobenius polynomial degree: %ld\n", fq_nmod_poly_degree(frobenius, ctx));
    
    fq_nmod_poly_gcd(splitting_poly, poly, frobenius, ctx);
    slong splitting_deg = fq_nmod_poly_degree(splitting_poly, ctx);
    
    //printf("GCD result degree: %ld (should be > 0 if there are roots in F_q)\n", splitting_deg);
    
    if (splitting_deg <= 0) {
        //printf("No roots found in F_q - this might be the problem!\n");
        goto cleanup;
    }
    
    if (splitting_deg == deg) {
        //printf("All roots are in F_q, proceeding with splitting\n");
        fq_nmod_poly_set(splitting_poly, poly, ctx);
    }
    
    if (splitting_deg == 1) {
        //printf("Found single linear factor\n");
        fq_nmod_t a, b, root;
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fq_nmod_init(root, ctx);
        
        fq_nmod_poly_get_coeff(a, splitting_poly, 1, ctx);
        fq_nmod_poly_get_coeff(b, splitting_poly, 0, ctx);
        
        if (!fq_nmod_is_zero(a, ctx)) {
            fq_nmod_inv(root, a, ctx);
            fq_nmod_neg(b, b, ctx);
            fq_nmod_mul(root, root, b, ctx);
            fq_nmod_roots_add(roots, root, 1, ctx);
            //printf("Extracted root: ");
            //fq_nmod_print_pretty(root, ctx);
            //printf("\n");
        }
        
        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(root, ctx);
    } else {
        //printf("Need random splitting for degree %ld polynomial\n", splitting_deg);
        
        flint_rand_t state;
        flint_rand_init(state);
        
        int success = 0;
        int max_attempts = 1000; 
        
        for (int attempt = 0; attempt < max_attempts && !success; attempt++) {
            //printf("Random splitting attempt %d\n", attempt + 1);
            
            fq_nmod_poly_t g, h, gcd_result;
            fq_nmod_poly_init(g, ctx);
            fq_nmod_poly_init(h, ctx);  
            fq_nmod_poly_init(gcd_result, ctx);
            
            fq_nmod_poly_randtest(g, state, splitting_deg, ctx);
            //printf("Random polynomial g degree: %ld\n", fq_nmod_poly_degree(g, ctx));
            
            fmpz_t exp;
            fmpz_init(exp);
            fmpz_sub_ui(exp, q, 1);
            fmpz_fdiv_q_2exp(exp, exp, 1);  // (q-1)/2
            
            //printf("Computing g^exp where exp = ");
            //fmpz_print(exp);
            //printf("\n");
            
            fq_nmod_poly_powmod_fmpz_binexp(h, g, exp, splitting_poly, ctx);
            
            // h = h - 1
            fq_nmod_t one;
            fq_nmod_init(one, ctx);
            fq_nmod_one(one, ctx);
            fq_nmod_poly_t one_poly;
            fq_nmod_poly_init(one_poly, ctx);
            fq_nmod_poly_set_fq_nmod(one_poly, one, ctx);
            
            fq_nmod_poly_sub(h, h, one_poly, ctx);
            
            //printf("h polynomial degree: %ld\n", fq_nmod_poly_degree(h, ctx));
            
            // gcd(h, splitting_poly)
            fq_nmod_poly_gcd(gcd_result, h, splitting_poly, ctx);
            slong gcd_deg = fq_nmod_poly_degree(gcd_result, ctx);
            
            //printf("GCD degree: %ld\n", gcd_deg);
            
            if (gcd_deg > 0 && gcd_deg < splitting_deg) {
                //printf("Successful split! Recursing...\n");
                
                our_fq_nmod_poly_roots(roots, gcd_result, with_multiplicity, ctx);
                
                fq_nmod_poly_t quotient;
                fq_nmod_poly_init(quotient, ctx);
                fq_nmod_poly_div(quotient, splitting_poly, gcd_result, ctx);
                our_fq_nmod_poly_roots(roots, quotient, with_multiplicity, ctx);
                fq_nmod_poly_clear(quotient, ctx);
                
                success = 1;
            } else {
                //printf("Split failed, trying again...\n");
            }
            
            fmpz_clear(exp);
            fq_nmod_clear(one, ctx);
            fq_nmod_poly_clear(one_poly, ctx);
            fq_nmod_poly_clear(g, ctx);
            fq_nmod_poly_clear(h, ctx);
            fq_nmod_poly_clear(gcd_result, ctx);
        }
        
        flint_rand_clear(state);
        
        if (!success) {
            //printf("Random splitting failed after %d attempts!\n", max_attempts);
        }
    }
    
cleanup:
    fmpz_clear(q);
    fq_nmod_poly_clear(x_to_q, ctx);
    fq_nmod_poly_clear(x, ctx);
    fq_nmod_poly_clear(frobenius, ctx);
    fq_nmod_poly_clear(splitting_poly, ctx);
    
    //printf("Total roots found: %ld\n", roots->num);
    return roots->num;
}

// Polynomial generation functions for testing
void generate_nmod_poly(nmod_poly_t poly, flint_rand_t state, slong degree, mp_limb_t p) {
    nmod_poly_zero(poly);
    mp_limb_t lead_coeff = n_randint(state, p - 1) + 1;
    nmod_poly_set_coeff_ui(poly, degree, lead_coeff);
    
    for (slong i = 0; i < degree; i++) {
        mp_limb_t coeff = n_randint(state, p);
        nmod_poly_set_coeff_ui(poly, i, coeff);
    }
}

void generate_fq_nmod_poly(fq_nmod_poly_t poly, flint_rand_t state, slong degree, const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_randtest_monic(poly, state, degree + 1, ctx);
}

// Benchmark functions

double benchmark_nmod_roots(slong degree, int num_tests) {
    printf("\n=== nmod_poly CZ Root Finding Test ===\n");
    
    nmod_poly_t poly;
    nmod_roots_t roots;
    flint_rand_t state;
    
    nmod_poly_init(poly, PRIME);
    nmod_roots_init(roots);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_nmod_poly(poly, state, degree, PRIME);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_nmod_poly_roots(roots, poly, 1);
        double end = get_time_roots();
        
        printf("nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    nmod_poly_clear(poly);
    nmod_roots_clear(roots);
    flint_rand_clear(state);
    
    return avg;
}

double benchmark_fq_nmod_roots(slong degree, slong extension, int num_tests) {
    printf("\n=== fq_nmod_poly CZ Root Finding Test (F_{%llu^%ld}) ===\n", 2, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t roots;
    flint_rand_t state;
    
    // Initialize finite field context
    fq_nmod_ctx_init_ui(ctx, 2, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(roots, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, ctx);
        double end = get_time_roots();
        
        printf("fq_nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("fq_nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

double benchmark_flint_fq_nmod_factor(slong degree, slong extension, int num_tests) {
    printf("\n=== FLINT fq_nmod_poly Factorization Test (F_{%llu^%ld}) ===\n", 2, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_poly_factor_t factors;
    flint_rand_t state;
    
    fq_nmod_ctx_init_ui(ctx, 2, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_poly_factor_init(factors, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        fq_nmod_poly_factor_clear(factors, ctx);
        fq_nmod_poly_factor_init(factors, ctx);
        
        fq_nmod_t lead_coeff;
        fq_nmod_init(lead_coeff, ctx);
        
        double start = get_time_roots();
        fq_nmod_poly_factor(factors, lead_coeff, poly, ctx);  // Correct parameter order
        double end = get_time_roots();
        
        slong linear_factors = 0;
        for (slong j = 0; j < factors->num; j++) {
            if (fq_nmod_poly_degree(factors->poly + j, ctx) == 1) {
                linear_factors += factors->exp[j];
            }
        }
        
        printf("FLINT fq_nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, linear_factors);
        total_time += (end - start);
        total_roots += linear_factors;
        
        fq_nmod_clear(lead_coeff, ctx);
    }
    
    double avg = total_time / num_tests;
    printf("FLINT fq_nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_poly_factor_clear(factors, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

// Test and verification functions

void test_fq_nmod_correctness() {
    printf("\n=== fq_nmod_poly Correctness Verification Test ===\n");
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t our_roots;
    flint_rand_t state;
    
    // Use F_5^2
    fq_nmod_ctx_init_ui(ctx, 5, 2, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(our_roots, ctx);
    flint_rand_init(state);
    
    printf("Test field: F_{5^2} = F_25\n");
    
    // Construct a simple polynomial: x^2 - 1 = (x-1)(x+1)
    fq_nmod_poly_zero(poly, ctx);
    fq_nmod_t one, neg_one;
    fq_nmod_init(one, ctx);
    fq_nmod_init(neg_one, ctx);
    fq_nmod_one(one, ctx);
    fq_nmod_set_ui(neg_one, 4, ctx);  // -1 in F_5 is 4
    
    fq_nmod_poly_set_coeff(poly, 2, one, ctx);      // x^2
    fq_nmod_poly_set_coeff(poly, 0, neg_one, ctx);  // -1 (in F_5, -1 = 4)
    
    printf("Test polynomial: x^2 - 1\n");
    printf("Expected roots: 1, 4 (in F_5)\n");
    
    our_roots->num = 0;
    slong num_roots = our_fq_nmod_poly_roots(our_roots, poly, 1, ctx);
    
    printf("Our found roots: ");
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(" ");
    }
    printf("(total %ld)\n", num_roots);
    
    // Verify root correctness
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_t value;
        fq_nmod_init(value, ctx);
        fq_nmod_poly_evaluate_fq_nmod(value, poly, our_roots->roots + i, ctx);
        printf("Verify f(");
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(") = ");
        fq_nmod_print_pretty(value, ctx);
        printf("\n");
        fq_nmod_clear(value, ctx);
    }
    
    fq_nmod_clear(one, ctx);
    fq_nmod_clear(neg_one, ctx);
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(our_roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
}

void run_unified_comparison() {
    printf("=== Fixed Version Unified CZ Root Finding Algorithm Test ===\n");
    printf("Comparing nmod_poly and fq_nmod_poly versions\n");
    printf("====================================\n");
    
    test_fq_nmod_correctness();
    
    slong degrees[] = {100, 500, 1000};
    int num_degrees = sizeof(degrees) / sizeof(degrees[0]);
    int num_tests = 3;
    
    printf("\n=== nmod_poly Test (F_%llu) ===\n", PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- Degree %ld ---\n", degrees[i]);
        benchmark_nmod_roots(degrees[i], num_tests);
    }
    
    printf("\n=== fq_nmod_poly Test (F_{%llu^2}) ===\n", SMALL_PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- Degree %ld ---\n", degrees[i]);
        double our_time = benchmark_fq_nmod_roots(degrees[i], 8, num_tests);
        double flint_time = benchmark_flint_fq_nmod_factor(degrees[i], 8, num_tests);
        
        double ratio = (flint_time > 0) ? our_time / flint_time : 0;
        printf("Ratio (ours/FLINT): %.2fx\n", ratio);
    }
}
