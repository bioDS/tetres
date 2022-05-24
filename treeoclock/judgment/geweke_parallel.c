#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


void compute_geweke_list(long n, long dm[][n], double geweke[n], double first_s, double first_e, double last)
{
    #pragma omp parallel for
        for(int i=10; i<n; i++)
        {
            long intersum = 0;
            long intersum_division = 0;
            long first_sum = 0;
            long last_sum = 0;

            for(int r=floor(i * first_s); r<floor(i * first_e); r++)
            {
                for(int s=floor(i * first_s); s<floor(i * first_e); s++)
                {
                    first_sum += dm[r][s] + dm[s][r];
//                    if(r < s)
//                    {
//                        first_sum += dm[r][s];
//                    }
//                    else
//                    {
//                        first_sum += dm[s][r];
//                    }
                }
            }
            for(int r=floor(i * (1 - last)); r<i; r++)
            {
                for(int s=floor(i * (1 - last)); s<i; s++)
                {
                    last_sum += dm[r][s] + dm[s][r];
//                    if(r < s)
//                    {
//                        last_sum += dm[r][s];
//                    }
//                    else
//                    {
//                        last_sum += dm[s][r];
//                    }
                }
            }
            for(int r=floor(i * first_s); r<floor(i * first_e); r++)
            {
                for(int s=floor(i * (1 - last)); s<i; s++)
                {
                        intersum += dm[r][s] + dm[s][r];
                        intersum_division++;
                }
            }

            double first_len = floor(i * first_e) - floor(i * first_s);
            double last_len = i - floor(i * (1 - last));
            first_len = pow(first_len, 2.0) - first_len;
            last_len = pow(last_len, 2.0) - last_len;
            if(first_len == 0){
                first_len = 1;
            }
            if(last_len == 0)
            {
                last_len = 1;
            }


//            (first_sum / first_len)
//            (last_sum / last_len)
//            (intersum / intersum_division)
            double result = 0;
            result = (double) fabs(((double) first_sum / (double) first_len) - ((double) last_sum / (double) last_len));
            // printf("%ld\n", i);


            result += (double) fabs(((double) first_sum / (double) first_len) - ((double) intersum / (double) intersum_division));
            result += (double) fabs(((double) last_sum / (double) last_len) - ((double) intersum / (double) intersum_division));

//            if(500 <= i && i <= 525)
//            {
//                   // printf("C: %f\n", sqrt(result));
//                   // printf("Result%d: %ld \t %f \t %ld \t %f \t %ld \t %ld\n", i, first_sum, first_len, last_sum, last_len,  intersum, intersum_division);
//                // printf("C: %f, %f\n", ((double) first_sum / (double) first_len), ((double) last_sum / (double) last_len));
//            }

            geweke[i] = sqrt(result);
//            floor(i * first_s)
//            floor(i * first_e)
//            floor(i * (1 - last))

        }

//        for(int i=0; i < n; i++)
//        {
//            for (int j=0; j < n; j++)
//            {
//                dm[i][j] = findpath_distance(&tree_set_one->trees[i], &tree_set_two->trees[j]);
//            }
//        }
}
