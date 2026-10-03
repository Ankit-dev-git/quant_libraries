#include <gtest/gtest.h>
#include "normal.h"
#include <cmath>

TEST(NormalTest, KnownValues){
    double pdf_0 = normal::norm_pdf(0.0);
    double pdf_1 = normal::norm_pdf(3, 3, 2);
    double cdf_0 = normal::norm_cdf(0);
    double cdf_1 = normal::norm_cdf(1.96);
    double cdf_2 = normal::norm_cdf(5, 3, 2);
    EXPECT_NEAR(pdf_0, 0.3989422804014327, 1e-14);
    EXPECT_NEAR(pdf_1, 0.19947114020071635, 1e-14);
    EXPECT_NEAR(cdf_0, 0.5, 1e-14);
    EXPECT_NEAR(cdf_1, 0.9750021048517795, 1e-14);
    EXPECT_NEAR(cdf_2, 0.8413447460685429, 1e-14);
}

TEST(NormalTest, Symmetry){
    double pdf_x = normal::norm_pdf(2);
    double pdf_mx = normal::norm_pdf(-2);
    double cdf_x = normal::norm_cdf(3);
    double cdf_mx = 1- normal::norm_cdf(-3);
    double cdf_x1 = normal::norm_cdf (3,1);
    double cdf_mx1 = 1 - normal::norm_cdf (-1,1);
    EXPECT_NEAR (pdf_x, pdf_mx, 1e-14);
    EXPECT_NEAR(cdf_x,  cdf_mx,  1e-14);
    EXPECT_NEAR(cdf_x1, cdf_mx1, 1e-14);
}

TEST(NormalTest, Derivative){
    double pdf_derv = (normal::norm_cdf(3+1e-5) - normal::norm_cdf(3-1e-5))/(2*1e-5);
    double pdf_dir = normal::norm_pdf(3);
    double pdf_derv_1 = (normal::norm_cdf(3+1e-5,1,2) - normal::norm_cdf(3-1e-5,1,2))/(2*1e-5);
    double pdf_dir_1 = normal::norm_pdf(3,1,2);
    EXPECT_NEAR(pdf_derv, pdf_dir, 1e-8);
    EXPECT_NEAR(pdf_derv_1, pdf_dir_1, 1e-8);
}

TEST(NormalTest, BadInput){
    EXPECT_THROW(normal::norm_pdf(1,0, -1), std::invalid_argument);
    EXPECT_THROW(normal::norm_pdf(1,0, 0), std::invalid_argument);
    EXPECT_THROW(normal::norm_cdf(1,0, -1), std::invalid_argument);
    EXPECT_THROW(normal::norm_cdf(1,0, 0), std::invalid_argument);
}

TEST(NormalTest, LeftTail){
    double cdf = normal::norm_cdf(-10);
    EXPECT_NEAR(cdf/7.619853024160593e-24, 1, 1e-12);
}