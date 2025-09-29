// dixon_gui.c - 主文件，不包含任何dixon头文件
// x86_64-w64-mingw32-gcc -O3 -march=native -fopenmp -static -o dixon_gui.exe *.c -mwindows -I../include -L../lib -lflint -lmpfr -lgmp -lm -lpthread -lstdc++ -lcomctl32 -lcomdlg32 -luser32 -lgdi32 -lshlwapi


// dixon_gui.c - 完全隔离版本，不包含任何FLINT头文件
#include <windows.h>
#include <libloaderapi.h>
#include <shlwapi.h>
#include <commctrl.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Window and control IDs
#define ID_TAB_CONTROL        1001
#define ID_BASIC_POLYS        1010
#define ID_BASIC_VARS         1011
#define ID_BASIC_FIELD        1012
#define ID_SOLVER_POLYS       1020
#define ID_SOLVER_FIELD       1021
#define ID_IDEAL_POLYS        1030
#define ID_IDEAL_VARS         1031
#define ID_IDEAL_GENERATORS   1032
#define ID_IDEAL_FIELD        1034
#define ID_RESULT_TEXT        1040
#define ID_BTN_COMPUTE        1050
#define ID_BTN_CLEAR          1051
#define ID_BTN_SAVE           1052
#define ID_BTN_LOAD           1053
#define ID_PROGRESS_BAR       1060
#define ID_STATUS_BAR         1061

// Custom messages
#define WM_COMPUTATION_DONE   (WM_USER + 1)
/*
// 包装函数声明 - 只使用基础C类型
extern int validate_field_size(const char *field_str, char *error_msg, int error_msg_size);
extern char* get_field_info(const char *field_str);
extern char* dixon_compute_basic(const char *polys_str, const char *vars_str, const char *field_str);
extern char* dixon_compute_solver(const char *polys_str, const char *field_str);
extern char* dixon_compute_ideal(const char *polys_str, const char *vars_str, 
                                const char *ideal_str, const char *field_str);
*/

// 动态加载的函数指针类型
typedef int (*validate_field_size_func)(const char *field_str, char *error_msg, int error_msg_size);
typedef char* (*get_field_info_func)(const char *field_str);
typedef char* (*dixon_compute_basic_func)(const char *polys_str, const char *vars_str, const char *field_str);
typedef char* (*dixon_compute_solver_func)(const char *polys_str, const char *field_str);
typedef char* (*dixon_compute_ideal_func)(const char *polys_str, const char *vars_str, 
                                         const char *ideal_str, const char *field_str);
// DLL相关全局变量

static validate_field_size_func validate_field_size = NULL;
static get_field_info_func get_field_info = NULL;
static dixon_compute_basic_func dixon_compute_basic = NULL;
static dixon_compute_solver_func dixon_compute_solver = NULL;
static dixon_compute_ideal_func dixon_compute_ideal = NULL;

static HMODULE g_dixonDll = NULL;
// 全局变量
HINSTANCE g_hInst;
HWND g_hMainWnd;
HWND g_hTabControl;
HWND g_hProgressBar;
HWND g_hStatusBar;
HWND g_hResultText;

// Controls for tabs
HWND g_hBasicPolys, g_hBasicVars, g_hBasicField;
HWND g_hSolverPolys, g_hSolverField;
HWND g_hIdealPolys, g_hIdealVars, g_hIdealGenerators, g_hIdealField;

// Buttons
HWND g_hBtnCompute, g_hBtnClear, g_hBtnSave, g_hBtnLoad;

// Threading
HANDLE g_hComputeThread = NULL;
BOOL g_bComputing = FALSE;

// 计算线程数据 - 只使用基础类型
typedef struct {
    char *polys_str;
    char *vars_str;
    char *field_str;
    char *ideal_str;
    int computation_mode;   // 0=basic dixon, 1=solver, 2=dixon with ideal
    char *result;
    BOOL success;
    double computation_time;
} ComputeData;

// 延迟加载DLL函数
static BOOL LoadDixonDLL() {
    if (g_dixonDll) return TRUE;  // 已经加载
    /*
    // 写入日志记录加载时机
    FILE *f = fopen("dll_load_log.txt", "w");
    if (f) {
        fprintf(f, "Attempting to load dixon_math.dll AFTER constructor\n");
        fclose(f);
    }
    */
    // 加载DLL
    g_dixonDll = LoadLibraryA("dixon.dll");
    if (!g_dixonDll) {
        char error[512];
        sprintf(error, "Failed to load dixon_math.dll\nError code: %lu\n\nMake sure:\n1. dixon_math.dll is in the same directory\n2. All FLINT DLLs are accessible", GetLastError());
        MessageBoxA(NULL, error, "DLL Load Error", MB_OK | MB_ICONERROR);
        return FALSE;
    }
	
    
    // 获取函数地址
    validate_field_size = (validate_field_size_func)GetProcAddress(g_dixonDll, "validate_field_size");
    get_field_info = (get_field_info_func)GetProcAddress(g_dixonDll, "get_field_info");
    dixon_compute_basic = (dixon_compute_basic_func)GetProcAddress(g_dixonDll, "dixon_compute_basic");
    dixon_compute_solver = (dixon_compute_solver_func)GetProcAddress(g_dixonDll, "dixon_compute_solver");
    dixon_compute_ideal = (dixon_compute_ideal_func)GetProcAddress(g_dixonDll, "dixon_compute_ideal");
    
    if (!validate_field_size || !get_field_info || !dixon_compute_basic || 
        !dixon_compute_solver || !dixon_compute_ideal) {
        MessageBoxA(NULL, "Failed to get function addresses from dixon_math.dll\n\nCheck if all functions are properly exported", "Function Import Error", MB_OK | MB_ICONERROR);
        FreeLibrary(g_dixonDll);
        g_dixonDll = NULL;
        return FALSE;
    }
    
    return TRUE;
}


// 获取控件文本
static char* GetEditText(HWND hEdit) {
    int len = GetWindowTextLengthA(hEdit);
    if (len == 0) return strdup("");
    
    char *buffer = malloc(len + 1);
    GetWindowTextA(hEdit, buffer, len + 1);
    return buffer;
}

// 设置控件文本
static void SetEditText(HWND hEdit, const char *text) {
    SetWindowTextA(hEdit, text ? text : "");
}

// 设置状态栏文本
static void SetStatusText(const char *text) {
    SendMessageA(g_hStatusBar, SB_SETTEXTA, 0, (LPARAM)text);
}

// 计算线程函数 - 只调用包装函数
static DWORD WINAPI ComputeThreadProc(LPVOID lpParam) {
    ComputeData *data = (ComputeData*)lpParam;
    clock_t start_time = clock();

	// 确保DLL已加载（这里才第一次触碰FLINT相关代码）
    if (!LoadDixonDLL()) {
        data->success = FALSE;
        PostMessage(g_hMainWnd, WM_COMPUTATION_DONE, 0, (LPARAM)data);
        return 0;
    }

    // 根据模式调用相应的包装函数
    if (data->computation_mode == 1) {
        // 多项式系统求解
        data->result = dixon_compute_solver(data->polys_str, data->field_str);
    } else if (data->computation_mode == 2) {
        // Dixon with ideal reduction
        data->result = dixon_compute_ideal(data->polys_str, data->vars_str, 
                                         data->ideal_str, data->field_str);
    } else {
        // 基础Dixon resultant
        data->result = dixon_compute_basic(data->polys_str, data->vars_str, data->field_str);
    }
    
    data->success = (data->result != NULL);
    
    // 计算时间
    clock_t end_time = clock();
    data->computation_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    // 通知主线程
    PostMessage(g_hMainWnd, WM_COMPUTATION_DONE, 0, (LPARAM)data);
    return 0;
}

// 统计逗号分隔项数量
static int count_comma_separated_items(const char *str) {
    if (!str || strlen(str) == 0) return 0;
    
    int count = 1;
    for (const char *p = str; *p; p++) {
        if (*p == ',') count++;
    }
    return count;
}

// 开始计算
static void StartComputation() {
    if (g_bComputing) return;
    
    // 获取当前选项卡
    int current_tab = TabCtrl_GetCurSel(g_hTabControl);
    
    ComputeData *data = malloc(sizeof(ComputeData));
    memset(data, 0, sizeof(ComputeData));
    
    if (current_tab == 0) {
        // 多项式系统求解器
        data->polys_str = GetEditText(g_hSolverPolys);
        data->vars_str = NULL;
        data->field_str = GetEditText(g_hSolverField);
        data->ideal_str = NULL;
        data->computation_mode = 1;
    } else if (current_tab == 1) {
        // 基础Dixon
        data->polys_str = GetEditText(g_hBasicPolys);
        data->vars_str = GetEditText(g_hBasicVars);
        data->field_str = GetEditText(g_hBasicField);
        data->ideal_str = NULL;
        data->computation_mode = 0;
    } else {
        // Dixon with ideal reduction
        data->polys_str = GetEditText(g_hIdealPolys);
        data->vars_str = GetEditText(g_hIdealVars);
        data->field_str = GetEditText(g_hIdealField);
        data->ideal_str = GetEditText(g_hIdealGenerators);
        data->computation_mode = 2;
    }
    
    // 验证基础输入
    if (strlen(data->polys_str) == 0 || strlen(data->field_str) == 0) {
        MessageBoxA(g_hMainWnd, "Please fill polynomials and field size!", "Input Error", MB_OK | MB_ICONERROR);
        goto cleanup_and_return;
    }
    
    // 使用包装函数验证字段大小
    char error_msg[512];
    if (!validate_field_size(data->field_str, error_msg, sizeof(error_msg))) {
        char full_error[1024];
        sprintf(full_error, "Invalid field size '%s'\n\n%s\n\nField size must be:\n"
                           "- A prime number (e.g., 257)\n"
                           "- A prime power (e.g., 256 = 2^8)\n"
                           "- In p^k format (e.g., 2^8, 3^5)", 
                           data->field_str, error_msg);
        MessageBoxA(g_hMainWnd, full_error, "Field Size Error", MB_OK | MB_ICONERROR);
        goto cleanup_and_return;
    }
    
    // 显示字段信息
    char *field_info = get_field_info(data->field_str);
    if (field_info) {
        SetStatusText(field_info);
        free(field_info);
    }
    
    // Dixon模式的额外验证
    if (data->computation_mode == 0 || data->computation_mode == 2) {
        if (!data->vars_str || strlen(data->vars_str) == 0) {
            MessageBoxA(g_hMainWnd, "Please specify variables to eliminate for Dixon resultant!", "Input Error", MB_OK | MB_ICONERROR);
            goto cleanup_and_return;
        }
        
        // 检查方程/变量数量
        int poly_count = count_comma_separated_items(data->polys_str);
        int var_count = count_comma_separated_items(data->vars_str);
        
        if (var_count != poly_count - 1) {
            char warning[512];
            sprintf(warning, "Warning: Dixon method requires eliminating exactly %d variables for %d equations!\r\n\r\nYou have %d variables to eliminate.\r\n\r\nContinue anyway?", 
                    poly_count - 1, poly_count, var_count);
            
            int result = MessageBoxA(g_hMainWnd, warning, "Variable Count Warning", MB_YESNO | MB_ICONWARNING);
            if (result != IDYES) {
                goto cleanup_and_return;
            }
        }
    }
    
    // ideal reduction模式的额外验证
    if (data->computation_mode == 2) {
        if (!data->ideal_str || strlen(data->ideal_str) == 0) {
            MessageBoxA(g_hMainWnd, "Please specify ideal generators for Dixon with ideal reduction!", "Input Error", MB_OK | MB_ICONERROR);
            goto cleanup_and_return;
        }
    }
    
    // 更新UI状态
    g_bComputing = TRUE;
    EnableWindow(g_hBtnCompute, FALSE);
    EnableWindow(g_hBtnSave, FALSE);
    ShowWindow(g_hProgressBar, SW_SHOW);
    SendMessage(g_hProgressBar, PBM_SETMARQUEE, TRUE, 50);
    SetStatusText("Computing...");
    
    // 启动计算线程
    DWORD threadId;
    g_hComputeThread = CreateThread(NULL, 0, ComputeThreadProc, data, 0, &threadId);
    
    if (!g_hComputeThread) {
        MessageBoxA(g_hMainWnd, "Failed to create computation thread!", "Error", MB_OK | MB_ICONERROR);
        g_bComputing = FALSE;
        EnableWindow(g_hBtnCompute, TRUE);
        ShowWindow(g_hProgressBar, SW_HIDE);
        SetStatusText("Ready");
        goto cleanup_and_return;
    }
    
    return;

cleanup_and_return:
    if (data->polys_str) free(data->polys_str);
    if (data->vars_str) free(data->vars_str);
    if (data->field_str) free(data->field_str);
    if (data->ideal_str) free(data->ideal_str);
    free(data);
}

// 处理计算完成
static void HandleComputationDone(ComputeData *data) {
    // 更新UI状态
    g_bComputing = FALSE;
    EnableWindow(g_hBtnCompute, TRUE);
    ShowWindow(g_hProgressBar, SW_HIDE);
    SendMessage(g_hProgressBar, PBM_SETMARQUEE, FALSE, 0);
    
    if (data->success && data->result) {
        // 构建显示文本
        size_t result_len = strlen(data->result);
        size_t msg_len = result_len + 512;
        char *display_text = malloc(msg_len);
        
        if (display_text) {
            const char *mode_names[] = {"Dixon Resultant", "Polynomial System Solver", "Dixon Resultant with Ideal Reduction"};
            snprintf(display_text, msg_len, 
                    "%s completed!\r\nTime: %.3f seconds\r\n\r\n%s", 
                    mode_names[data->computation_mode], data->computation_time, data->result);
            
            SetEditText(g_hResultText, display_text);
            free(display_text);
        } else {
            SetEditText(g_hResultText, data->result);
        }
        
        SetStatusText("Computation completed");
        EnableWindow(g_hBtnSave, TRUE);
    } else {
        // 处理失败情况
        const char *error_msgs[] = {
            "Dixon resultant computation failed! Please check your input.\r\n\r\nNote: Number of variables to eliminate must equal (number of equations - 1).",
            "Polynomial system solving failed! Please check your input.\r\n\r\nNote: Number of equations must equal number of variables.",
            "Dixon resultant with ideal reduction failed! Please check your input.\r\n\r\nNote: Check polynomial syntax and ideal generators."
        };
        
        SetEditText(g_hResultText, error_msgs[data->computation_mode]);
        SetStatusText("Computation failed");
    }
    
    // 清理数据
    if (data->polys_str) free(data->polys_str);
    if (data->vars_str) free(data->vars_str);
    if (data->field_str) free(data->field_str);
    if (data->ideal_str) free(data->ideal_str);
    if (data->result) free(data->result);
    free(data);
    
    // 关闭线程句柄
    if (g_hComputeThread) {
        CloseHandle(g_hComputeThread);
        g_hComputeThread = NULL;
    }
}

// 清除结果
static void ClearResults() {
    SetEditText(g_hResultText, "");
    SetStatusText("Ready");
    EnableWindow(g_hBtnSave, FALSE);
}

// 保存结果
static void SaveResults() {
    OPENFILENAMEA ofn;
    char szFile[260] = "dixon_result.txt";
    
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_hMainWnd;
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "Text Files\0*.txt\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
    
    if (GetSaveFileNameA(&ofn)) {
        char *text = GetEditText(g_hResultText);
        
        FILE *file = fopen(szFile, "w");
        if (file) {
            fprintf(file, "%s", text);
            fclose(file);
            SetStatusText("Results saved");
        } else {
            MessageBoxA(g_hMainWnd, "Failed to save file!", "Error", MB_OK | MB_ICONERROR);
        }
        
        free(text);
    }
}

// 加载文件
static void LoadFile() {
    OPENFILENAMEA ofn;
    char szFile[260] = "";
    
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_hMainWnd;
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "Data Files\0*.dat\0Text Files\0*.txt\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
    
    if (GetOpenFileNameA(&ofn)) {
        FILE *file = fopen(szFile, "r");
        if (file) {
            char line[4096];
            int data_line_count = 0;
            char *field_str = NULL, *polys_str = NULL, *vars_str = NULL, *ideal_str = NULL;
            
            while (fgets(line, sizeof(line), file)) {
                // 去除换行符
                line[strcspn(line, "\n\r")] = '\0';
                
                // 跳过空行和注释
                if (strlen(line) == 0 || line[0] == '#') {
                    continue;
                }
                
                // 处理数据行
                if (data_line_count == 0) {
                    field_str = strdup(line);
                } else if (data_line_count == 1) {
                    polys_str = strdup(line);
                } else if (data_line_count == 2) {
                    vars_str = strdup(line);
                } else if (data_line_count == 3) {
                    ideal_str = strdup(line);
                    break;
                }
                
                data_line_count++;
                if (data_line_count >= 4) break;
            }
            
            fclose(file);
            
            // 根据当前选项卡填充字段
            int current_tab = TabCtrl_GetCurSel(g_hTabControl);
            
            if (current_tab == 0) {
                // 多项式系统求解器选项卡
                if (field_str) SetEditText(g_hSolverField, field_str);
                if (polys_str) SetEditText(g_hSolverPolys, polys_str);
            } else if (current_tab == 1) {
                // 基础Dixon结果选项卡
                if (field_str) SetEditText(g_hBasicField, field_str);
                if (polys_str) SetEditText(g_hBasicPolys, polys_str);
                if (vars_str) SetEditText(g_hBasicVars, vars_str);
            } else if (current_tab == 2) {
                // Dixon with Ideal Reduction选项卡
                if (field_str) SetEditText(g_hIdealField, field_str);
                if (polys_str) SetEditText(g_hIdealPolys, polys_str);
                if (vars_str) SetEditText(g_hIdealVars, vars_str);
                if (ideal_str) SetEditText(g_hIdealGenerators, ideal_str);
            }
            
            // 显示成功消息
            char success_msg[512];
            snprintf(success_msg, sizeof(success_msg), 
                    "File loaded successfully!\n\nLoaded %d data lines:\n%s%s%s%s", 
                    data_line_count,
                    field_str ? "- Field size\n" : "",
                    polys_str ? "- Polynomials\n" : "",
                    vars_str ? "- Variables\n" : "",
                    ideal_str ? "- Ideal generators\n" : "");
            MessageBoxA(g_hMainWnd, success_msg, "Load Success", MB_OK | MB_ICONINFORMATION);
            
            SetStatusText("File loaded successfully");
            
            // 清理分配的字符串
            if (field_str) free(field_str);
            if (polys_str) free(polys_str);
            if (vars_str) free(vars_str);
            if (ideal_str) free(ideal_str);
            
        } else {
            MessageBoxA(g_hMainWnd, "Cannot open file!", "Error", MB_OK | MB_ICONERROR);
        }
    }
}

// 创建选项卡控件
static void CreateTabControls(HWND hParent) {
    int y = 10;
    int edit_height = 22;
    int label_height = 16;
    int margin = 8;
    int client_width = 560;
    
    // 根据选项卡索引创建不同控件
    int tab_index = (int)GetWindowLongPtrA(hParent, GWLP_USERDATA);
    
    if (tab_index == 0) {
        // 多项式系统求解器选项卡
        CreateWindowA("STATIC", "Polynomial System (comma separated):", WS_CHILD | WS_VISIBLE,
                     10, y, 320, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hSolverPolys = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                      10, y, client_width, edit_height * 2, hParent, (HMENU)ID_SOLVER_POLYS, g_hInst, NULL);
        y += edit_height * 2 + margin;
        
        CreateWindowA("STATIC", "Field size (examples: 257, 256, 2^8, 3^5):", WS_CHILD | WS_VISIBLE,
                     10, y, 280, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hSolverField = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                      10, y, 120, edit_height, hParent, (HMENU)ID_SOLVER_FIELD, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Note: Number of equations must equal number of variables (n*n)", 
                     WS_CHILD | WS_VISIBLE, 10, y, 480, label_height, hParent, NULL, g_hInst, NULL);
        
        SendMessageA(g_hSolverPolys, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: x + y - 3, 2*x - y  or  x^2 + y^2 - 1, x + y - 1");
        SendMessageA(g_hSolverField, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: 257, 2^8, 256");
        
    } else if (tab_index == 1) {
        // 基础Dixon选项卡
        CreateWindowA("STATIC", "Polynomials (comma separated):", WS_CHILD | WS_VISIBLE,
                     10, y, 280, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hBasicPolys = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                     10, y, client_width, edit_height, hParent, (HMENU)ID_BASIC_POLYS, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Variables to eliminate (comma separated):", WS_CHILD | WS_VISIBLE,
                     10, y, 320, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hBasicVars = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                    10, y, 280, edit_height, hParent, (HMENU)ID_BASIC_VARS, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Field size (examples: 257, 256, 2^8, 3^5):", WS_CHILD | WS_VISIBLE,
                     10, y, 280, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hBasicField = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                     10, y, 120, edit_height, hParent, (HMENU)ID_BASIC_FIELD, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Note: Number of variables to eliminate must equal (number of equations - 1)", 
                     WS_CHILD | WS_VISIBLE, 10, y, 500, label_height, hParent, NULL, g_hInst, NULL);
        
        SendMessageA(g_hBasicPolys, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: x + y + z, x*y + y*z + z*x, x*y*z + 1");
        SendMessageA(g_hBasicVars, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: x, y");
        SendMessageA(g_hBasicField, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: 257, 2^8, 256");
        
    } else {
        // Dixon with ideal reduction选项卡
        CreateWindowA("STATIC", "Polynomials (comma separated):", WS_CHILD | WS_VISIBLE,
                     10, y, 280, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hIdealPolys = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                     10, y, client_width, edit_height, hParent, (HMENU)ID_IDEAL_POLYS, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Variables to eliminate:", WS_CHILD | WS_VISIBLE,
                     10, y, 180, label_height, hParent, NULL, g_hInst, NULL);
        CreateWindowA("STATIC", "Field size:", WS_CHILD | WS_VISIBLE,
                     300, y, 80, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hIdealVars = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                    10, y, 280, edit_height, hParent, (HMENU)ID_IDEAL_VARS, g_hInst, NULL);
        g_hIdealField = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                     300, y, 120, edit_height, hParent, (HMENU)ID_IDEAL_FIELD, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Ideal generators (comma separated):", WS_CHILD | WS_VISIBLE,
                     10, y, 300, label_height, hParent, NULL, g_hInst, NULL);
        y += label_height + 2;
        g_hIdealGenerators = CreateWindowA("EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                          10, y, client_width, edit_height, hParent, (HMENU)ID_IDEAL_GENERATORS, g_hInst, NULL);
        y += edit_height + margin;
        
        CreateWindowA("STATIC", "Note: Number of variables to eliminate must equal (number of equations - 1)", 
                     WS_CHILD | WS_VISIBLE, 10, y, 500, label_height, hParent, NULL, g_hInst, NULL);
        
        SendMessageA(g_hIdealPolys, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5");
        SendMessageA(g_hIdealVars, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: a4");
        SendMessageA(g_hIdealGenerators, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3");
        SendMessageA(g_hIdealField, EM_SETCUEBANNER, TRUE, (LPARAM)L"e.g.: 257, 2^8, 256");
    }
}

// 窗口过程
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_CREATE: {
        // 初始化通用控件
        INITCOMMONCONTROLSEX icex;
        icex.dwSize = sizeof(INITCOMMONCONTROLSEX);
        icex.dwICC = ICC_TAB_CLASSES | ICC_PROGRESS_CLASS | ICC_BAR_CLASSES;
        InitCommonControlsEx(&icex);
        
        // 创建选项卡控件
        g_hTabControl = CreateWindowA(WC_TABCONTROLA, "",
                                     WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS,
                                     10, 10, 580, 200, hWnd, (HMENU)ID_TAB_CONTROL, g_hInst, NULL);
        
        // 添加选项卡
        TCITEMA tie;
        tie.mask = TCIF_TEXT;
        tie.pszText = "Polynomial System Solver";
        TabCtrl_InsertItem(g_hTabControl, 0, &tie);
        
        tie.pszText = "Dixon Resultant";
        TabCtrl_InsertItem(g_hTabControl, 1, &tie);
        
        tie.pszText = "Dixon with Ideal Reduction";
        TabCtrl_InsertItem(g_hTabControl, 2, &tie);
        
        // 创建选项卡内容窗口
        RECT rcTab;
        GetClientRect(g_hTabControl, &rcTab);
        TabCtrl_AdjustRect(g_hTabControl, FALSE, &rcTab);
        
        HWND hSolverTab = CreateWindowA("STATIC", "", WS_CHILD | WS_VISIBLE,
                                       rcTab.left, rcTab.top, rcTab.right - rcTab.left, rcTab.bottom - rcTab.top,
                                       g_hTabControl, NULL, g_hInst, NULL);
        
        HWND hBasicTab = CreateWindowA("STATIC", "", WS_CHILD,
                                      rcTab.left, rcTab.top, rcTab.right - rcTab.left, rcTab.bottom - rcTab.top,
                                      g_hTabControl, NULL, g_hInst, NULL);
        
        HWND hIdealTab = CreateWindowA("STATIC", "", WS_CHILD,
                                      rcTab.left, rcTab.top, rcTab.right - rcTab.left, rcTab.bottom - rcTab.top,
                                      g_hTabControl, NULL, g_hInst, NULL);
        
        // 存储选项卡句柄和索引
        SetWindowLongPtrA(hSolverTab, GWLP_USERDATA, 0);
        SetWindowLongPtrA(hBasicTab, GWLP_USERDATA, 1);
        SetWindowLongPtrA(hIdealTab, GWLP_USERDATA, 2);
        
        CreateTabControls(hSolverTab);
        CreateTabControls(hBasicTab);
        CreateTabControls(hIdealTab);
        
        // 创建按钮
        g_hBtnCompute = CreateWindowA("BUTTON", "Compute", WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
                                     10, 220, 80, 28, hWnd, (HMENU)ID_BTN_COMPUTE, g_hInst, NULL);
        
        g_hBtnClear = CreateWindowA("BUTTON", "Clear", WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
                                   100, 220, 80, 28, hWnd, (HMENU)ID_BTN_CLEAR, g_hInst, NULL);
        
        g_hBtnSave = CreateWindowA("BUTTON", "Save", WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
                                  190, 220, 80, 28, hWnd, (HMENU)ID_BTN_SAVE, g_hInst, NULL);
        
        g_hBtnLoad = CreateWindowA("BUTTON", "Load", WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
                                  280, 220, 80, 28, hWnd, (HMENU)ID_BTN_LOAD, g_hInst, NULL);
        
        EnableWindow(g_hBtnSave, FALSE);
        
        // 创建进度条
        g_hProgressBar = CreateWindowA(PROGRESS_CLASSA, NULL,
                                      WS_CHILD | PBS_MARQUEE,
                                      10, 258, 580, 18, hWnd, (HMENU)ID_PROGRESS_BAR, g_hInst, NULL);
        
        // 创建结果文本区域
        CreateWindowA("STATIC", "Computation Results:", WS_CHILD | WS_VISIBLE,
                     10, 285, 150, 16, hWnd, NULL, g_hInst, NULL);
        
        g_hResultText = CreateWindowA("EDIT", "",
                                     WS_CHILD | WS_VISIBLE | WS_BORDER | WS_VSCROLL | WS_HSCROLL |
                                     ES_MULTILINE | ES_READONLY,
                                     10, 305, 580, 140, hWnd, (HMENU)ID_RESULT_TEXT, g_hInst, NULL);
        
        // 设置等宽字体
        HFONT hFont = CreateFontA(12, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, ANSI_CHARSET,
                                 OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY,
                                 FIXED_PITCH | FF_MODERN, "Consolas");
        SendMessage(g_hResultText, WM_SETFONT, (WPARAM)hFont, TRUE);
        
        // 创建状态栏
        g_hStatusBar = CreateWindowA(STATUSCLASSNAMEA, NULL,
                                    WS_CHILD | WS_VISIBLE | SBARS_SIZEGRIP,
                                    0, 0, 0, 0, hWnd, (HMENU)ID_STATUS_BAR, g_hInst, NULL);
        
        SetStatusText("Ready - Complete FLINT Isolation (No FLINT headers in GUI)");
        
        break;
    }
    
    case WM_COMMAND: {
        switch (LOWORD(wParam)) {
        case ID_BTN_COMPUTE:
            StartComputation();
            break;
        case ID_BTN_CLEAR:
            ClearResults();
            break;
        case ID_BTN_SAVE:
            SaveResults();
            break;
        case ID_BTN_LOAD:
            LoadFile();
            break;
        }
        break;
    }
    
    case WM_NOTIFY: {
        LPNMHDR pnmh = (LPNMHDR)lParam;
        if (pnmh->idFrom == ID_TAB_CONTROL && pnmh->code == TCN_SELCHANGE) {
            // 处理选项卡选择变化
            int sel = TabCtrl_GetCurSel(g_hTabControl);
            
            // 隐藏所有选项卡内容窗口
            HWND hChild = GetWindow(g_hTabControl, GW_CHILD);
            while (hChild) {
                ShowWindow(hChild, SW_HIDE);
                hChild = GetWindow(hChild, GW_HWNDNEXT);
            }
            
            // 显示选中的选项卡内容
            hChild = GetWindow(g_hTabControl, GW_CHILD);
            int index = 0;
            while (hChild && index <= sel) {
                if (index == sel) {
                    ShowWindow(hChild, SW_SHOW);
                    break;
                }
                hChild = GetWindow(hChild, GW_HWNDNEXT);
                index++;
            }
        }
        break;
    }
    
    case WM_COMPUTATION_DONE: {
        HandleComputationDone((ComputeData*)lParam);
        break;
    }
    
    case WM_SIZE: {
        // 重新调整状态栏大小
        SendMessage(g_hStatusBar, WM_SIZE, wParam, lParam);
        break;
    }
    
    case WM_CLOSE: {
        if (g_bComputing) {
            int result = MessageBoxA(hWnd, "Computation is in progress. Are you sure you want to exit?", "Confirm Exit", 
                                    MB_YESNO | MB_ICONQUESTION);
            if (result != IDYES) {
                return 0;
            }
        }
        
        if (g_hComputeThread) {
            TerminateThread(g_hComputeThread, 0);
            CloseHandle(g_hComputeThread);
        }
        
        DestroyWindow(hWnd);
        break;
    }
    
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    
    return 0;
}

// 入口点
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
    wchar_t exePath[MAX_PATH];
    GetModuleFileNameW(NULL, exePath, MAX_PATH);
    PathRemoveFileSpecW(exePath);
    
    wchar_t dllPath[MAX_PATH];
    PathCombineW(dllPath, exePath, L"dll");
    
    if (GetFileAttributesW(dllPath) != INVALID_FILE_ATTRIBUTES) {
        SetDefaultDllDirectories(LOAD_LIBRARY_SEARCH_DEFAULT_DIRS | LOAD_LIBRARY_SEARCH_USER_DIRS);
        AddDllDirectory(dllPath);
    }
	LoadDixonDLL();
	
    g_hInst = hInstance;
    
    // 注册窗口类
    WNDCLASSEXA wcex;
    wcex.cbSize = sizeof(WNDCLASSEXA);
    wcex.style = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc = WndProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(hInstance, IDI_APPLICATION);
    wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName = NULL;
    wcex.lpszClassName = "DixonGUIWindow";
    wcex.hIconSm = LoadIcon(wcex.hInstance, IDI_APPLICATION);
    
    if (!RegisterClassExA(&wcex)) {
        MessageBoxA(NULL, "Window class registration failed!", "Error", MB_OK | MB_ICONERROR);
        return 1;
    }
    
    // 创建主窗口
    g_hMainWnd = CreateWindowA("DixonGUIWindow", "Dixon Resultant & Polynomial System Solver - Complete FLINT Isolation",
                              WS_OVERLAPPEDWINDOW,
                              CW_USEDEFAULT, 0, 620, 500, NULL, NULL, hInstance, NULL);
    
    if (!g_hMainWnd) {
        MessageBoxA(NULL, "Window creation failed!", "Error", MB_OK | MB_ICONERROR);
        return 1;
    }
    
    ShowWindow(g_hMainWnd, nCmdShow);
    UpdateWindow(g_hMainWnd);
    
    // 消息循环
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    
    return (int)msg.wParam;
}