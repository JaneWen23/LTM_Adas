/**
 * 
 * Function：LTM实现
 * Author:joker.mao
 * 
*/

#ifndef _LTM_H_
#define _LTM_H_

#include <assert.h>
#include <math.h>
#include <iostream>

#define MAX_DIV_NUMS     256
#define MAX_H_V_DIV_NUMS 16
#define ALPHA_A          0.0

#define HIST_DOMAIN    1  // 0 is for linear, 1 is for log
#define HIST_CEILING   1  // 0: disabled; 1: enabled
#define LUMA_MIN       0.1f // min should be greater than 0.0f
#define LUMA_MAX       1.0f // should be greater than LUMA_MIN
#define CEILING_K      5.0f

#define SAFE_FREE_ARR(mem) if (mem) {delete[] mem;}

float int255_to_ln(int val){
    // scale 255
    float result = logf(((float)val/255.0f)*(LUMA_MAX-LUMA_MIN) + LUMA_MIN); // debug
    return result;
}

float normalize_ln(float lnValF){
    const float lnMin = logf(LUMA_MIN);
    const float lnMax = logf(LUMA_MAX);
    return (lnValF - lnMin) / (lnMax - lnMin);
}

float denormalize_ln(float lnValFNorm){
    const float lnMin = logf(LUMA_MIN);
    const float lnMax = logf(LUMA_MAX);
    return lnValFNorm * (lnMax - lnMin) + lnMin;
}

int exp_ln_255(float lnValF){
    return (int)((expf(lnValF) - LUMA_MIN) / (LUMA_MAX-LUMA_MIN) * 255.0f);
}

void test_int_and_ln(){
    for (int i = 0; i < 256; ++i){
        float lnValF = int255_to_ln(i);
        float lnValFNorm = normalize_ln(lnValF);
        printf("%d, ln=%f, normLn = %f, exp(ln) = %d\n", i, lnValF, lnValFNorm, exp_ln_255(lnValF));
    }
}

template<typename _TI, typename _TO>
class Ltm {
private:
    float *div_img_historm_bins_[MAX_H_V_DIV_NUMS][MAX_H_V_DIV_NUMS];
    _TO   *div_img_map_lut_[MAX_H_V_DIV_NUMS][MAX_H_V_DIV_NUMS];
    _TI   *div_img_data_[MAX_H_V_DIV_NUMS][MAX_H_V_DIV_NUMS];
    int   ave_bins_[MAX_H_V_DIV_NUMS][MAX_H_V_DIV_NUMS];
    int   sd_bins_[MAX_H_V_DIV_NUMS][MAX_H_V_DIV_NUMS];
    int   max_sd_ = 0;
    _TO   *out_;

    int width_  = 0;
    int height_ = 0;
    int bin_nums_ = 0;
    int bin_pixel_nums_ = 0;
    int hor_div_nums_;
    int ver_div_nums_;
    int hor_pixel_nums_per_bin_;
    int ver_pixel_nums_per_bin_;

    int in_max_;
    int out_max_;
    float phase_s = 3.5;
    float phase_d = 0;
public:
    Ltm(const _TI * const src, _TO* out, int w, int h, int h_bins, int v_bins, int input_max, int output_max) {
        out_ = out;

        assert(w % h_bins == 0 && h % v_bins == 0 && src != nullptr && out != nullptr);

        width_  = w;
        height_ = h;
        hor_div_nums_ = h_bins;
        ver_div_nums_ = v_bins;
        hor_pixel_nums_per_bin_ = w / h_bins;
        ver_pixel_nums_per_bin_ = h / v_bins;
        in_max_ = input_max;
        out_max_ = output_max;

        bin_nums_ = hor_div_nums_ * ver_div_nums_;
        bin_pixel_nums_ = hor_pixel_nums_per_bin_ * ver_pixel_nums_per_bin_;
        phase_d = sqrtf((height_ * height_) + (width_ * width_));

        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                //malloc data
                div_img_historm_bins_[bin_idx][bin_idy] = new float[input_max];
                div_img_map_lut_[bin_idx][bin_idy]      = new _TO[input_max];
                div_img_data_[bin_idx][bin_idy]         = new _TI[bin_pixel_nums_];
                assert(div_img_historm_bins_[bin_idx][bin_idy]);
                assert(div_img_map_lut_[bin_idx][bin_idy]);
                assert(div_img_data_[bin_idx][bin_idy]);
                memset(div_img_historm_bins_[bin_idx][bin_idy], 0, sizeof(float) * input_max);
                _TI val;
                int sum = 0;
                ave_bins_[bin_idx][bin_idy] = 0;
                sd_bins_[bin_idx][bin_idy] = 0;

                int max = 0;
                int min = in_max_;
                //copy div image
                for (int idy = 0; idy < ver_pixel_nums_per_bin_; ++idy) {
                    for (int idx = 0; idx < hor_pixel_nums_per_bin_; ++idx) {
                        val = div_img_data_[bin_idx][bin_idy][idy * hor_pixel_nums_per_bin_ + idx] = \
                            src[(bin_idy * ver_pixel_nums_per_bin_ + idy) * w + (bin_idx * hor_pixel_nums_per_bin_ + idx)];
                        // val is in linear domain
                        //=== val in log domain:========
                        #if HIST_DOMAIN == 1
                        float tmp = normalize_ln(int255_to_ln(val)) * 255.0f;
                        val = (_TI)tmp;                   
                        #endif
                        //===============================
                        ++div_img_historm_bins_[bin_idx][bin_idy][val];
                        sum += val;
                        if (min > val) min = val;
                        if (max < val) max = val;
                    }
                }
                // apply hist ceiling ====
                #if HIST_CEILING == 1
                float k = CEILING_K;
                HistCeiling(k, bin_idx, bin_idy);
                #endif
                //========================

                sd_bins_[bin_idx][bin_idy] = max - min;
                if (sd_bins_[bin_idx][bin_idy] > max_sd_) max_sd_ = sd_bins_[bin_idx][bin_idy];
                ave_bins_[bin_idx][bin_idy] = sum / bin_pixel_nums_;
                //init end
            }
        }
    };

    ~Ltm() {
        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                //malloc data
                SAFE_FREE_ARR(div_img_historm_bins_[bin_idx][bin_idy]); 
                SAFE_FREE_ARR(div_img_map_lut_[bin_idx][bin_idy]);
                SAFE_FREE_ARR(div_img_data_[bin_idx][bin_idy]);
            }
        }
    };

    void HistormBanlance() {
        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                float alpha = ALPHA_A * (1 - exp(-(max_sd_ - sd_bins_[bin_idx][bin_idy])));
                int numOfSmpls = 0;
                #if HIST_CEILING == 1
                    for(int i = 0; i < 256; ++i){
                        numOfSmpls += div_img_historm_bins_[bin_idx][bin_idy][i];
                    }
                #else
                    numOfSmpls = bin_pixel_nums_;
                #endif
                div_img_historm_bins_[bin_idx][bin_idy][0] /= numOfSmpls;
                div_img_map_lut_[bin_idx][bin_idy][0] = alpha * (div_img_historm_bins_[bin_idx][bin_idy][0]  * out_max_) \
                                + (1 - alpha) * (0 * out_max_ / in_max_);
                for (int i = 1; i < in_max_; ++i) {
                    div_img_historm_bins_[bin_idx][bin_idy][i] = div_img_historm_bins_[bin_idx][bin_idy][i] / numOfSmpls\
                                                             + div_img_historm_bins_[bin_idx][bin_idy][i - 1];
                    

                    // original:
                    // div_img_map_lut_[bin_idx][bin_idy][i] = alpha * (div_img_historm_bins_[bin_idx][bin_idy][i]  * out_max_) \
                    //             + (1 - alpha) * (i * out_max_ / in_max_);

                    //======= log domain ===========
                    #if HIST_DOMAIN == 1
                    float tmp = normalize_ln(int255_to_ln(i)) * 255.0f; // debug
                    int linearHist = (int)tmp;
                    alpha = 1.0;
                    #else
                    int linearHist = i;
                    alpha = 1.0;
                    #endif
                    //==============================
                    div_img_map_lut_[bin_idx][bin_idy][i] = alpha * (div_img_historm_bins_[bin_idx][bin_idy][i] * 255) \
                                + (1 - alpha) * (linearHist * out_max_ / in_max_);


                    // if (bin_idy == 0 && bin_idx == 0) {
                    //    //printf("%d, %f,   ", i, div_img_historm_bins_[bin_idx][bin_idy][i] * 255);
                    //    printf("%d   %d\n", (int)div_img_map_lut_[bin_idx][bin_idy][i], static_cast<int>(exp_ln_255(denormalize_ln(div_img_map_lut_[bin_idx][bin_idy][i] / 255.0f))));
                    // }
                }
            }
        }
    }

    void HistCeiling(float k, int bin_idx, int bin_idy){
        int avgSmplPerBin = bin_pixel_nums_ / in_max_; // in_max_ is total bin number
        float tolerance = bin_pixel_nums_ * 0.05;
        int ceiling = (int)(avgSmplPerBin * k); // n = total number of samples / total bin number = avg samples_per_bin
        int trimmings = (int)tolerance + 1;
        int totalSmpls = 0;
        while(trimmings > tolerance){
            trimmings = 0;
            totalSmpls = 0;
            
            for (int i = 0; i < in_max_; ++i){
                if((int)div_img_historm_bins_[bin_idx][bin_idy][i] > ceiling){
                    trimmings += (int)div_img_historm_bins_[bin_idx][bin_idy][i] - ceiling;
                    div_img_historm_bins_[bin_idx][bin_idy][i] = ceiling;
                }
                totalSmpls += (int)div_img_historm_bins_[bin_idx][bin_idy][i];
            }
            ceiling = (int)(k * totalSmpls / in_max_);
            // debug:
            // printf("block(%d,%d), ceiling = %d, trimmings = %d\n", bin_idy, bin_idx, ceiling, trimmings);
        }
        // debug:
            printf("block(%d,%d), ceiling = %d, trimmings = %d\n", bin_idy, bin_idx, ceiling, trimmings);
        // for (int i = 0; i < 256; ++i){
        //     printf("hist[%d] = %d\t", i, (int)div_img_historm_bins_[bin_idx][bin_idy][i]);
        // }
    }

    int GetPixelMapVal(int val, int ave, int x, int y) {
        double ws_wd_sum = 0;
        double ws_wd_map_sum = 0;

        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                //===== log domain ========
                #if HIST_DOMAIN == 1
                float valForLut = normalize_ln(int255_to_ln(val)) * 255.0f; // debug
                #else
                float valForLut = val;
                #endif
                //=========================
                float ws = 1; //exp(-(fabs(valForLut - ave_bins_[bin_idx][bin_idy]) / in_max_) / phase_s); // debug: val to lnVal
                int hc = bin_idx * hor_pixel_nums_per_bin_ + (hor_pixel_nums_per_bin_ >> 1);
                int vc = bin_idy * ver_pixel_nums_per_bin_ + (ver_pixel_nums_per_bin_ >> 1);
                float wd = 1; //exp(-(sqrt((x - hc) * (x - hc) + (y - vc) *(y - vc))) / phase_d);
                ws_wd_sum += (ws * wd);
                //ws_wd_map_sum += (div_img_map_lut_[bin_idx][bin_idy][(int)lnValF] * ws * wd); // equation (12) of Duan2010; val to lnVal
                ws_wd_map_sum += (div_img_map_lut_[bin_idx][bin_idy][(int)valForLut] * ws * wd); // equation (12) of Duan2010;
            }
        }
        
        //===== log domain =========
        #if HIST_DOMAIN == 1
        int rt = static_cast<int>(exp_ln_255(denormalize_ln(ws_wd_map_sum/ (255.0f*ws_wd_sum)))); // debug, exp(Luma)
        #else
        int rt = static_cast<int>(ws_wd_map_sum / ws_wd_sum);
        #endif
        //==========================
        if (rt >= out_max_) {
            rt = out_max_ - 1;
        }
        return rt;
    }

    void RunLtmMap() {
        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                for (int idy = 0; idy < ver_pixel_nums_per_bin_; ++idy) {
                    for (int idx = 0; idx < hor_pixel_nums_per_bin_; ++idx) {
                        auto& val = div_img_data_[bin_idx][bin_idy][idy * hor_pixel_nums_per_bin_ + idx];
                        //==== local only =========
                        //val = div_img_map_lut_[bin_idx][bin_idy][val];
                        //=========================
                        val = GetPixelMapVal(val, ave_bins_[bin_idx][bin_idy], (bin_idx * hor_pixel_nums_per_bin_ + idx), (bin_idy * ver_pixel_nums_per_bin_ + idy));
                    }
                }
            }
        }
    }

    void CombineImg() {
        for (int bin_idy = 0; bin_idy < ver_div_nums_; ++bin_idy) {
            for (int bin_idx = 0; bin_idx < hor_div_nums_; ++bin_idx) {
                //copy div image
                for (int idy = 0; idy < ver_pixel_nums_per_bin_; ++idy) {
                    for (int idx = 0; idx < hor_pixel_nums_per_bin_; ++idx) {
                        out_[(bin_idy * ver_pixel_nums_per_bin_ + idy) * width_ + (bin_idx * hor_pixel_nums_per_bin_ + idx)] \
                            = div_img_data_[bin_idx][bin_idy][idy * hor_pixel_nums_per_bin_ + idx];
                    }
                }
            }
        }
    }

    void Run() {
        HistormBanlance();
        RunLtmMap();
        CombineImg();
    };

    const float* GetDivImgHis(int idx, int idy) const {
        assert((idx < hor_div_nums_) && (idy < ver_div_nums_));
        return div_img_historm_bins_[idx][idy];
    }

    _TI* GetDivImgPtr(int idx, int idy) {
        assert((idx < hor_div_nums_) && (idy < ver_div_nums_));
        return div_img_data_[idx][idy];
    }  
    _TO* GetLutMapPtr(int idx, int idy) {
        assert((idx < hor_div_nums_) && (idy < ver_div_nums_));
        return div_img_map_lut_[idx][idy];
    }
};



#endif