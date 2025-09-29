// dixon_wrapper.c - 完全隔离的包装模块，包含所有FLINT相关代码
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 跨平台DLL/共享库导出定义  __declspec(dllexport)
#ifdef DLL_EXPORT
	#define EXPORT
#else
	#define EXPORT
#endif

// 包含所有FLINT头文件 - 在DLL路径设置之后安全加载
#include <flint/fmpz_factor.h>
#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "polynomial_system_solver.h"

// 内部质数幂检查函数
static int check_prime_power_internal(const fmpz_t n, fmpz_t prime, unsigned long *power) {
    if (fmpz_cmp_ui(n, 1) <= 0) return 0;
    
    // 检查是否直接是质数
    if (fmpz_is_probabprime(n)) {
        fmpz_set(prime, n);
        *power = 1;
        return 1;
    }
    
    // 使用FLINT的因式分解功能
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    
    // 检查是否只有一个质因子（即质数的幂）
    if (factors->num == 1) {
        fmpz_set(prime, factors->p + 0);  // 第一个（也是唯一的）质因子
        *power = factors->exp[0];         // 对应的幂
        fmpz_factor_clear(factors);
        return 1;
    }
    
    fmpz_factor_clear(factors);
    return 0;
}

// 解析字段大小的内部函数
static int parse_field_size_internal(const char *field_str, fmpz_t prime, unsigned long *power) {
    if (!field_str || strlen(field_str) == 0) {
        return 0;
    }
    
    // 检查是否包含'^'（幂记号）
    const char *caret = strchr(field_str, '^');
    if (caret) {
        // 格式：p^k
        char *prime_str = malloc(caret - field_str + 1);
        strncpy(prime_str, field_str, caret - field_str);
        prime_str[caret - field_str] = '\0';
        
        // 解析质数部分
        fmpz_t p;
        fmpz_init(p);
        int success = fmpz_set_str(p, prime_str, 10);
        if (success != 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // 解析幂部分
        char *endptr;
        unsigned long k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // 检查p是否为质数
        if (!fmpz_is_probabprime(p)) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        fmpz_set(prime, p);
        *power = k;
        fmpz_clear(p);
        free(prime_str);
        return 1;
    } else {
        // 格式：直接数字（可能是p或p^k）
        fmpz_t field_size;
        fmpz_init(field_size);
        int success = fmpz_set_str(field_size, field_str, 10);
        if (success != 0) {
            fmpz_clear(field_size);
            return 0;
        }
        
        // 检查是否为质数幂
        int result = check_prime_power_internal(field_size, prime, power);
        fmpz_clear(field_size);
        return result;
    }
}

// 公共接口：验证字段大小
EXPORT int validate_field_size(const char *field_str, char *error_msg, int error_msg_size) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    int result = parse_field_size_internal(field_str, prime, &power);
    
    if (!result) {
        if (error_msg && error_msg_size > 0) {
            snprintf(error_msg, error_msg_size, "Failed to parse '%s' as a valid field size", field_str);
        }
    } else {
        // 检查质数是否太大
        if (!fmpz_fits_si(prime)) {
            if (error_msg && error_msg_size > 0) {
                snprintf(error_msg, error_msg_size, "Prime is too large (must fit in machine word)");
            }
            result = 0;
        }
    }
    
    fmpz_clear(prime);
    return result;
}

// 公共接口：获取字段信息
EXPORT char* get_field_info(const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Invalid field");
    }
    
    char *info = malloc(256);
    if (!info) {
        fmpz_clear(prime);
        return NULL;
    }
    
    if (power == 1) {
        // 质数域
        char *prime_str = fmpz_get_str(NULL, 10, prime);
        snprintf(info, 256, "Field: F_%s (prime field)", prime_str);
        //free(prime_str);
    } else {
        // 扩域
        unsigned long prime_ui = fmpz_get_ui(prime);
        unsigned long field_size = 1;
        for (unsigned long i = 0; i < power; i++) {
            field_size *= prime_ui;
        }
        snprintf(info, 256, "Field: F_%lu^%lu (extension field, size %lu, generator 't')", 
                prime_ui, power, field_size);
    }
    
    fmpz_clear(prime);
    return info;
}

// 公共接口：基础Dixon resultant计算
EXPORT char* dixon_compute_basic(const char *polys_str, const char *vars_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // 初始化有限域上下文
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // 调用Dixon resultant函数
    char *result = dixon_str(polys_str, vars_str, ctx);
    
    // 清理
    //fq_nmod_ctx_clear(ctx);
    //fmpz_clear(prime);
    
    //if (!result) {
        //return strdup("Error: Dixon resultant computation failed");
    //}
    
    return result;
}


// Simple solution formatting in computation thread - GCC COMPATIBLE VERSION WITH WINDOWS LINE ENDINGS
static char* format_solutions_simple(const polynomial_solutions_t *sols, const fq_nmod_ctx_t ctx) {
    if (!sols) {
        return strdup("Solution structure is null\r\n");
    }
    
    char *buffer = malloc(16384);  // 16KB buffer
    if (!buffer) {
        return strdup("Memory allocation failed\r\n");
    }
    
    strcpy(buffer, "\r\n=== POLYNOMIAL SYSTEM SOLUTIONS ===\r\n");
    
    if (!sols->is_valid) {
        strcat(buffer, "Solving failed");
        if (sols->error_message) {
            strcat(buffer, ": ");
            strncat(buffer, sols->error_message, 500);
        }
        strcat(buffer, "\r\n");
        return buffer;
    }
    
    if (sols->has_no_solutions == 1) {
        strcat(buffer, "System has no solutions over the finite field\r\n");
        return buffer;
    }
    
    if (sols->has_no_solutions == -1) {
        strcat(buffer, "Solving failed: polynomial system dimension greater than zero, please use Dixon resultant elimination\r\n");
        return buffer;
    }
    
    if (sols->num_variables == 0) {
        strcat(buffer, "No variables\r\n");
        return buffer;
    }
    
    if (sols->num_solution_sets == 0) {
        strcat(buffer, "No solutions found\r\n");
        return buffer;
    }
    
    char temp[512];
    sprintf(temp, "Found %ld complete solution set(s):\r\n", sols->num_solution_sets);
    strcat(buffer, temp);
    
    // Format solutions with simple error handling
    for (slong set = 0; set < sols->num_solution_sets && set < 5; set++) {
        sprintf(temp, "\r\nSolution set %ld:\r\n", set + 1);
        strcat(buffer, temp);
        
        for (slong var = 0; var < sols->num_variables && var < 10; var++) {
            if (!sols->variable_names || !sols->variable_names[var]) {
                continue;
            }
            
            sprintf(temp, "  %s = ", sols->variable_names[var]);
            strcat(buffer, temp);
            
            if (!sols->solutions_per_var) {
                strcat(buffer, "ERROR: solutions_per_var is null\r\n");
                continue;
            }
            
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            
            if (num_sols == 0) {
                strcat(buffer, "no solution\r\n");
            } else if (num_sols == 1) {
                if (!sols->solution_sets || !sols->solution_sets[set] || 
                    !sols->solution_sets[set][var] || !sols->solution_sets[set][var][0]) {
                    strcat(buffer, "ERROR: Invalid solution pointer\r\n");
                } else {
                    // Simple approach: try conversion, use fallback if it fails
                    char *sol_str = NULL;
                    
                    // Try the pretty print conversion
                    sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][0], ctx);
                    
                    // Check if conversion was successful
                    if (sol_str && strlen(sol_str) > 0 && strlen(sol_str) < 500) {
                        strncat(buffer, sol_str, 200);
                        strcat(buffer, "\r\n");
                    } else {
                        // Fallback: For prime fields, try to show as simple integer
                        if (fq_nmod_ctx_degree(ctx) == 1) {
                            // For prime fields, extract the coefficient
                            fmpz_t coeff;
                            fmpz_init(coeff);
                            fq_nmod_get_fmpz(coeff, sols->solution_sets[set][var][0], ctx);
                            
                            if (fmpz_fits_si(coeff)) {
                                slong val = fmpz_get_si(coeff);
                                sprintf(temp, "%ld\r\n", val);
                                strcat(buffer, temp);
                            } else {
                                strcat(buffer, "[large integer value]\r\n");
                            }
                            fmpz_clear(coeff);
                        } else {
                            strcat(buffer, "[field extension element]\r\n");
                        }
                    }
                    
                    // Clean up sol_str if it was allocated
                    if (sol_str) {
                        ;//free(sol_str); // This cause a segment fault
                    }
                }
            } else if (num_sols > 1 && num_sols <= 5) {
                strcat(buffer, "{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) strcat(buffer, ", ");
                    
                    if (sols->solution_sets && sols->solution_sets[set] && 
                        sols->solution_sets[set][var] && sols->solution_sets[set][var][sol]) {
                        
                        // For multiple solutions, use simple fallback
                        if (fq_nmod_ctx_degree(ctx) == 1) {
                            fmpz_t coeff;
                            fmpz_init(coeff);
                            fq_nmod_get_fmpz(coeff, sols->solution_sets[set][var][sol], ctx);
                            
                            if (fmpz_fits_si(coeff)) {
                                slong val = fmpz_get_si(coeff);
                                sprintf(temp, "%ld", val);
                                strcat(buffer, temp);
                            } else {
                                strcat(buffer, "[large_val]");
                            }
                            fmpz_clear(coeff);
                        } else {
                            strcat(buffer, "[ext_field]");
                        }
                    } else {
                        strcat(buffer, "NULL");
                    }
                }
                strcat(buffer, "}\r\n");
            } else {
                sprintf(temp, "%ld solutions (too many to display)\r\n", num_sols);
                strcat(buffer, temp);
            }
        }
    }
    
    strcat(buffer, "\r\n=== Solution Complete ===\r\n");
    return buffer;
}

// 公共接口：多项式系统求解
EXPORT char* dixon_compute_solver(const char *polys_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // 初始化有限域上下文
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // 调用多项式系统求解器
    polynomial_solutions_t *solutions = solve_polynomial_system_string(polys_str, ctx);
    char *result = format_solutions_simple(solutions, ctx);
    // 清理
    //fq_nmod_ctx_clear(ctx);
    //fmpz_clear(prime);
    
    //if (!result) {
        //return strdup("Error: Polynomial system solving failed");
    //}
    
    return result;
}

// 公共接口：Dixon with ideal reduction计算
EXPORT char* dixon_compute_ideal(const char *polys_str, const char *vars_str, 
                         const char *ideal_str, const char *field_str) {
    fmpz_t prime;
    unsigned long power;
    fmpz_init(prime);
    
    if (!parse_field_size_internal(field_str, prime, &power)) {
        fmpz_clear(prime);
        return strdup("Error: Invalid field size");
    }
    
    // 初始化有限域上下文
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, power, "t");
    
    // 调用Dixon with ideal reduction函数
    char *result = dixon_with_ideal_reduction_str(polys_str, vars_str, ideal_str, ctx);
    
    // 清理
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(prime);
    
    if (!result) {
        return strdup("Error: Dixon resultant with ideal reduction failed");
    }
    
    return result;
}


