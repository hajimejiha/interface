#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx2,tune=native")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

#define PI 3.14159265359

double rho_D[200][200][200]; //グローバル変数
 
int main(int argc, char *argv[]){
    //変更するパラメータ
    constexpr double dx = 1.0; //メッシュ間隔
    constexpr double K = 10.0; //粗視化パラメータ

    constexpr int grid_x = 250; //グリッド数
    constexpr int grid_y = 250;
    constexpr int grid_z = 250;

    constexpr double x_center = 100.0; //中心位置
    constexpr double y_center = 100.0;

    constexpr double z_0 = 50.0; //グリッド最下点
    
    const char* const dump_name_1 = "mydump_"; //lammpsダンプファイル名
    const char* const dump_name_2 = ".dmp";
    
    const char* const rho_name_1 = "rho_D_"; //出力ファイル名
    const char* const rho_name_2 = ".txt";
    
    const char* const dump_data = "%*d %d %*lf %lf %lf %lf %*lf %*lf %*lf %*lf"; //lammpsダンプファイルのフォーマット

    //コンパイル時に確定する定数
    constexpr double x_0 = x_center - (grid_x/2.0)*dx;
    constexpr double y_0 = y_center - (grid_y/2.0)*dx;

    constexpr double x_1 = x_0 - dx * K;
    constexpr double x_2 = x_center + (grid_x/2)*dx + dx * K;
    constexpr double y_1 = y_0 - dx * K;
    constexpr double y_2 = y_center + (grid_y/2)*dx + dx * K;
    constexpr double z_1 = z_0 - dx * K;
    constexpr double z_2 = z_0 + grid_z*dx + dx * K;

    constexpr double C1 = K * dx;
    constexpr double C2 = C1 / PI;    
    constexpr double C3 = pow(2.0*C1, -3.0);
    
    //プログラム実行時に確定する定数
    const int start_step = atoi(argv[1]);
    const int dump_step = atoi(argv[2]);
    const int file_num = atoi(argv[3]);  
    const int parallel_num = atoi(argv[4]);
    
    //変数
    int i, j, k, step, line_num, a_type;
    double X, Y, Z, temp_x, temp_y, temp_z, length;
    int count = 0;
    
    char dump_name[16];
    char rho_name[18];
    
    char buf[100];
    constexpr int buf_size = sizeof(buf);
    
    FILE *fp;
    
    for(int num = 0; num < file_num; ++num){
        clock_t start = clock();
        
        //0で初期化
        rho_D[grid_x-1][grid_y-1][grid_z-1] = {0.0};
        
        step = start_step + num * dump_step;
        printf("\nNow Step:%d", step);
        fflush(stdout);
        
        //ダンプファイル読み込み
        sprintf(dump_name, "%s%d%s", dump_name_1, step, dump_name_2);        
        fp = fopen(dump_name, "r");
        
        //開かなかった場合
        if(fp == NULL){
            printf("\nError: File can not be oppend\n");
            return 1;
        }
        
        line_num = 0;
        
        while(fgets(buf, buf_size, fp)){
            line_num++;
            
            //ダンプファイルの座標データまでcontinue
            if(line_num < 10) continue;
            
            sscanf(buf, dump_data, &a_type, &X, &Y, &Z);
                        
            if(a_type == 1){
                if(x_2 > X && X > x_1 && y_2 > Y && Y > y_1 && z_2 > Z && Z > z_1){
                    count++;
                    #pragma omp parallel for private(length, j, k, temp_x, temp_y, temp_z) num_threads(parallel_num)
                    for(i = 0; i < grid_x; ++i){
                        temp_x = i*dx+x_0-X;
                        for(j = 0; j < grid_y; ++j){
                            temp_y = j*dx+y_0-Y;
                            for(k = 0; k < grid_z; ++k){
                                temp_z = k*dx+z_0-Z;
                                
                                length = sqrt((temp_x * temp_x)
                                            +(temp_y * temp_y)
                                            +(temp_z * temp_z));                                
                                            
                                if(length < C1){
                                    rho_D[i][j][k] += C3 * (1.0+cos(C2*temp_x)) * (1.0+cos(C2*temp_y)) * (1.0+cos(C2*temp_z));
                                }
                            }
                        }
                    }
                }
            }
            //printf("\nLine:%d", line_num); //debug
            //fflush(stdout);
        }
        
        fclose(fp);
        
        //　密度データ書き込み
        sprintf(rho_name, "%s%d%s", rho_name_1, step, rho_name_2);        
        fp = fopen(rho_name, "w");
        
        //　開かなかった場合
        if(fp == NULL){
            printf("Error: File can not be oppend\n");
            return 1;
        }
        
        fprintf(fp, "step:%d atoms:%d dx:%lf K:%lf", step, count, dx, K);
        
        for(i = 0; i < grid_x; ++i){
            for(j = 0; j < grid_y; ++j){
                for(k = 0; k < grid_z; ++k){
                    fprintf(fp, "\n%6.3lf %6.3lf %6.3lf %10.8lf", i*dx+x_0, j*dx+y_0, k*dx+z_0, rho_D[i][j][k]);
                }
            }
        }
        
        fclose(fp);
        
        count = 0;
        
        clock_t end = clock();
        const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        printf("\ntime: %lf[sec]\n", time);
        fflush(stdout);
    }
    return 0;
}
