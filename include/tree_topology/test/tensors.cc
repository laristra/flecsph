#include "gtest/gtest.h"

#include <iomanip>
#include <iostream>
#include <log.h>

#include "tensor.h"

namespace flecsi {
namespace execution {
void
driver(int, char **) {}
} // namespace execution
} // namespace flecsi

TEST(tensors, sanity) {

  using namespace std;
  using namespace flecsi;
  using namespace tensor_indices;

  cout << "--- Generic tensor: ---" << endl;
  using gen_tensor_t = tensor_u<double, symmetry_type::generic, 3, 3>;
  gen_tensor_t G{0};
  G[xy] = 1.2;
  G[yz] = 2.3;
  G[xz] = 1.3;
  G(1, 1) = 2.0;
  G[xx] = 0.1;
  cout << "size of G: " << G.size() << endl;
  cout << "tensor G: " << endl << G << endl;

  cout << "--- Symmetric tensor: ---" << endl;
  using sym_tensor_t = tensor_u<double, symmetry_type::symmetric, 3, 3>;
  sym_tensor_t S{0};
  S[xy] = 1.2;
  S[yz] = 2.3;
  S[xz] = 1.3;
  S(1, 1) = 2.0;
  S[xx] = 0.1;
  cout << "size of S: " << S.size() << endl;
  cout << "tensor S: " << endl << S << endl;
  cout << "accessing via the [xy.. etc.] operator: " << endl
       << "S[xx]  S[xy]  S[xz]     ||" << setw(3) << S[xx] << " " << setw(3)
       << S[xy] << " " << setw(3) << S[xz] << " "
       << "||" << endl
       << "S[yx]  S[yy]  S[yz]  =  ||" << setw(3) << S[yx] << " " << setw(3)
       << S[yy] << " " << setw(3) << S[yz] << " "
       << "||" << endl
       << "S[zx]  S[zy]  S[zz]     ||" << setw(3) << S[zx] << " " << setw(3)
       << S[zy] << " " << setw(3) << S[zz] << " "
       << "||" << endl;
  cout << "accessing via the '(i,j)' operator:" << endl;
  for(int i = 0; i < 3; ++i) {
    cout << "S(" << i << ",0) S(" << i << ",1) S(" << i << ",2)";
    if(i == 1)
      cout << "  =  ||";
    else
      cout << "     ||";
    for(int j = 0; j < 3; ++j)
      cout << setw(3) << S(i, j) << " ";
    cout << "||" << endl;
  }

  cout << "--- Generic tensor, rank 3: ---" << endl;
  using gen2_tensor3_t = tensor_u<double, symmetry_type::generic, 2, 2, 2>;
  gen2_tensor3_t G3{0};
  G3[xxx] = 111;
  G3[xxy] = 112;
  G3[xyx] = 121;
  G3[xyy] = 122;
  G3(1, 0, 0) = 211;
  G3(1, 0, 1) = 212;
  G3(1, 1, 0) = 221;
  G3(1, 1, 1) = 222;
  cout << "size of G3: " << G3.size() << endl;
  cout << "tensor G3: " << endl << G3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: " << endl;
  cout << "G3[xxx] = " << G3[xxx] << "; G3[xxy] = " << G3[xxy] << endl;
  cout << "G3[xyx] = " << G3[xyx] << "; G3[xyy] = " << G3[xyy] << endl;
  cout << "G3[yxx] = " << G3[yxx] << "; G3[yxy] = " << G3[yxy] << endl;
  cout << "G3[yyx] = " << G3[yyx] << "; G3[yyy] = " << G3[yyy] << endl;

  cout << "--- Generic tensor, rank 3, dimension 3: ---" << endl;
  using gen3_tensor3_t = tensor_u<double, symmetry_type::generic, 3, 3, 3>;
  gen3_tensor3_t Q3{0};
  Q3[xxx] = 111;
  Q3[xxy] = 112;
  Q3[xyx] = 121;
  Q3[xyy] = 122;
  Q3(1, 0, 0) = 211;
  Q3(1, 0, 1) = 212;
  Q3(1, 1, 0) = 221;
  Q3(1, 1, 1) = 222;
  cout << "size of Q3: " << Q3.size() << endl;
  cout << "tensor Q3: " << endl << Q3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: " << endl;
  cout << "Q3[xxx] = " << Q3[xxx] << "; Q3[xxy] = " << Q3[xxy] << endl;
  cout << "Q3[xyx] = " << Q3[xyx] << "; Q3[xyy] = " << Q3[xyy] << endl;
  cout << "Q3[yxx] = " << Q3[yxx] << "; Q3[yxy] = " << Q3[yxy] << endl;
  cout << "Q3[yyx] = " << Q3[yyx] << "; Q3[yyy] = " << Q3[yyy] << endl;

  cout << "--- Symmetric tensor, rank 3: ---" << endl;
  using sym2_tensor3_t = tensor_u<double, symmetry_type::symmetric, 2, 2, 2>;
  sym2_tensor3_t S3{0};
  S3[xxx] = 111;
  S3[xxy] = 112;
  S3[xyy] = 122;
  S3(1, 1, 1) = 222;
  cout << "size of S3: " << S3.size() << endl;
  cout << "tensor S3: " << endl << S3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: " << endl;
  cout << "S3[xxx] = " << S3[xxx] << "; S3[xxy] = " << S3[xxy] << endl;
  cout << "S3[xyx] = " << S3[xyx] << "; S3[xyy] = " << S3[xyy] << endl;
  cout << "S3[yxx] = " << S3[yxx] << "; S3[yxy] = " << S3[yxy] << endl;
  cout << "S3[yyx] = " << S3[yyx] << "; S3[yyy] = " << S3[yyy] << endl;

  cout << "--- Symmetric tensor, rank 3, dimension 3: ---" << endl;
  using sym3_tensor3_t = tensor_u<double, symmetry_type::symmetric, 3, 3, 3>;
  sym3_tensor3_t Z3{0};
  Z3[xxx] = 111;
  Z3(_y, _y, _y) = 222;
  Z3(2, 2, 2) = 333;
  Z3[xxy] = 112;
  Z3[xyy] = 122;
  Z3[xxz] = 113;
  Z3[xzz] = 133;
  Z3[xyz] = 123;
  Z3[yyz] = 223;
  Z3[yzz] = 233;
  cout << "size of Z3: " << Z3.size() << endl;
  cout << "tensor Z3: " << endl << Z3 << endl;
  cout << "accessing via the [xyz.. etc.] operator: " << endl;
  cout << "Z3[xxx] = " << Z3[xxx] << "; Z3[xxy] = " << Z3[xxy]
       << "; Z3[xxz] = " << Z3[xxz] << endl;
  cout << "Z3[xyx] = " << Z3[xyx] << "; Z3[xyy] = " << Z3[xyy]
       << "; Z3[xyz] = " << Z3[xyz] << endl;
  cout << "Z3[xzx] = " << Z3[xzx] << "; Z3[xzy] = " << Z3[xzy]
       << "; Z3[xzz] = " << Z3[xzz] << endl;
  cout << endl;
  cout << "Z3[yxx] = " << Z3[yxx] << "; Z3[yxy] = " << Z3[yxy]
       << "; Z3[yxz] = " << Z3[yxz] << endl;
  cout << "Z3[yyx] = " << Z3[yyx] << "; Z3[yyy] = " << Z3[yyy]
       << "; Z3[yyz] = " << Z3[yyz] << endl;
  cout << "Z3[yzx] = " << Z3[yzx] << "; Z3[yzy] = " << Z3[yzy]
       << "; Z3[yzz] = " << Z3[yzz] << endl;
  cout << endl;
  cout << "Z3[zxx] = " << Z3[zxx] << "; Z3[zxy] = " << Z3[zxy]
       << "; Z3[zxz] = " << Z3[zxz] << endl;
  cout << "Z3[zyx] = " << Z3[zyx] << "; Z3[zyy] = " << Z3[zyy]
       << "; Z3[zyz] = " << Z3[zyz] << endl;
  cout << "Z3[zzx] = " << Z3[zzx] << "; Z3[zzy] = " << Z3[zzy]
       << "; Z3[zzz] = " << Z3[zzz] << endl;

  cout << "--- Generic tensor, rank 4, dimension 2: ---" << endl;
  using gen2_tensor4_t = tensor_u<double, symmetry_type::generic, 2, 2, 2, 2>;
  gen2_tensor4_t Q4{0};
  Q4[xxxx] = 1111;
  Q4[xxxy] = 1112;
  Q4[xxyx] = 1121;
  Q4[xxyy] = 1122;
  Q4(0, 1, 0, 0) = 1211;
  Q4(0, 1, 0, 1) = 1212;
  Q4(0, 1, 1, 0) = 1221;
  Q4(0, 1, 1, 1) = 1222;
  Q4(1, 0, 0, 0) = 2111;
  Q4(1, 0, 0, 1) = 2112;
  Q4(1, 0, 1, 0) = 2121;
  Q4(1, 0, 1, 1) = 2122;
  Q4(1, 1, 0, 0) = 2211;
  Q4(1, 1, 0, 1) = 2212;
  Q4(1, 1, 1, 0) = 2221;
  Q4(1, 1, 1, 1) = 2222;
  cout << "size of Q4: " << Q4.size() << endl;
  cout << "tensor Q4: " << endl << Q4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: " << endl;
  cout << "Q4[xxxx] = " << Q4[xxxx] << "; Q4[xxxy] = " << Q4[xxxy] << "; "
       << "Q4[xxyx] = " << Q4[xxyx] << "; Q4[xxyy] = " << Q4[xxyy] << endl;
  cout << "Q4[xyxx] = " << Q4[xyxx] << "; Q4[xyxy] = " << Q4[xyxy] << "; "
       << "Q4[xyyx] = " << Q4[xyyx] << "; Q4[xyyy] = " << Q4[xyyy] << endl;
  cout << "Q4[yxxx] = " << Q4[yxxx] << "; Q4[yxxy] = " << Q4[yxxy] << "; "
       << "Q4[yxyx] = " << Q4[yxyx] << "; Q4[yxyy] = " << Q4[yxyy] << endl;
  cout << "Q4[yyxx] = " << Q4[yyxx] << "; Q4[yyxy] = " << Q4[yyxy] << "; "
       << "Q4[yyyx] = " << Q4[yyyx] << "; Q4[yyyy] = " << Q4[yyyy] << endl;

  cout << "--- Generic tensor, rank 4, dimension 3: ---" << endl;
  using gen3_tensor4_t = tensor_u<double, symmetry_type::generic, 3, 3, 3, 3>;
  gen3_tensor4_t R4{0};
  R4[xxxx] = 1111;
  R4[xxxy] = 1112;
  R4[xxxz] = 1113;
  R4[xxyx] = 1121;
  R4[xxyy] = 1122;
  R4[xxyz] = 1123;
  R4[xxzx] = 1131;
  R4[xxzy] = 1132;
  R4[xxzz] = 1133;

  R4[xyxx] = 1211;
  R4[xyxy] = 1212;
  R4[xyxz] = 1213;
  R4[xyyx] = 1221;
  R4[xyyy] = 1222;
  R4[xyyz] = 1223;
  R4[xyzx] = 1231;
  R4[xyzy] = 1232;
  R4[xyzz] = 1233;

  R4[xzxx] = 1311;
  R4[xzxy] = 1312;
  R4[xzxz] = 1313;
  R4[xzyx] = 1321;
  R4[xzyy] = 1322;
  R4[xzyz] = 1323;
  R4[xzzx] = 1331;
  R4[xzzy] = 1332;
  R4[xzzz] = 1333;

  R4(1, 0, 0, 0) = 2111;
  R4(1, 0, 0, 1) = 2112;
  R4(1, 0, 0, 2) = 2113;
  R4(1, 0, 1, 0) = 2121;
  R4(1, 0, 1, 1) = 2122;
  R4(1, 0, 1, 2) = 2123;
  R4(1, 0, 2, 0) = 2131;
  R4(1, 0, 2, 1) = 2132;
  R4(1, 0, 2, 2) = 2133;

  R4[yyxx] = 2211;
  R4[yyxy] = 2212;
  R4[yyxz] = 2213;
  R4[yyyx] = 2221;
  R4[yyyy] = 2222;
  R4[yyyz] = 2223;
  R4[yyzx] = 2231;
  R4[yyzy] = 2232;
  R4[yyzz] = 2233;

  R4[yzxx] = 2311;
  R4[yzxy] = 2312;
  R4[yzxz] = 2313;
  R4[yzyx] = 2321;
  R4[yzyy] = 2322;
  R4[yzyz] = 2323;
  R4[yzzx] = 2331;
  R4[yzzy] = 2332;
  R4[yzzz] = 2333;

  R4[zxxx] = 3111;
  R4[zxxy] = 3112;
  R4[zxxz] = 3113;
  R4[zxyx] = 3121;
  R4[zxyy] = 3122;
  R4[zxyz] = 3123;
  R4[zxzx] = 3131;
  R4[zxzy] = 3132;
  R4[zxzz] = 3133;

  R4[zyxx] = 3211;
  R4[zyxy] = 3212;
  R4[zyxz] = 3213;
  R4[zyyx] = 3221;
  R4[zyyy] = 3222;
  R4[zyyz] = 3223;
  R4[zyzx] = 3231;
  R4[zyzy] = 3232;
  R4[zyzz] = 3233;

  R4[zzxx] = 3311;
  R4[zzxy] = 3312;
  R4[zzxz] = 3313;
  R4[zzyx] = 3321;
  R4[zzyy] = 3322;
  R4[zzyz] = 3323;
  R4[zzzx] = 3331;
  R4[zzzy] = 3332;
  R4[zzzz] = 3333;

  cout << "size of R4: " << R4.size() << endl;
  cout << "tensor R4: " << endl << R4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: " << endl;
  cout << "R4[xxxx] = " << R4[xxxx] << "; "
       << "R4[xxxy] = " << R4[xxxy] << "; "
       << "R4[xxxz] = " << R4[xxxz] << endl
       << "R4[xxyx] = " << R4[xxyx] << "; "
       << "R4[xxyy] = " << R4[xxyy] << "; "
       << "R4[xxyz] = " << R4[xxyz] << endl
       << "R4[xxzx] = " << R4[xxzx] << "; "
       << "R4[xxzy] = " << R4[xxzy] << "; "
       << "R4[xxzz] = " << R4[xxzz] << endl
       << endl;
  cout << "R4[xyxx] = " << R4[xyxx] << "; "
       << "R4[xyxy] = " << R4[xyxy] << "; "
       << "R4[xyxz] = " << R4[xyxz] << endl
       << "R4[xyyx] = " << R4[xyyx] << "; "
       << "R4[xyyy] = " << R4[xyyy] << "; "
       << "R4[xyyz] = " << R4[xyyz] << endl
       << "R4[xyzx] = " << R4[xyzx] << "; "
       << "R4[xyzy] = " << R4[xyzy] << "; "
       << "R4[xyzz] = " << R4[xyzz] << endl
       << endl;
  cout << "R4[xzxx] = " << R4[xzxx] << "; "
       << "R4[xzxy] = " << R4[xzxy] << "; "
       << "R4[xzxz] = " << R4[xzxz] << endl
       << "R4[xzyx] = " << R4[xzyx] << "; "
       << "R4[xzyy] = " << R4[xzyy] << "; "
       << "R4[xzyz] = " << R4[xzyz] << endl
       << "R4[xzzx] = " << R4[xzzx] << "; "
       << "R4[xzzy] = " << R4[xzzy] << "; "
       << "R4[xzzz] = " << R4[xzzz] << endl
       << endl;
  cout << endl;

  cout << "R4[yxxx] = " << R4[yxxx] << "; "
       << "R4[yxxy] = " << R4[yxxy] << "; "
       << "R4[yxxz] = " << R4[yxxz] << endl
       << "R4[yxyx] = " << R4[yxyx] << "; "
       << "R4[yxyy] = " << R4[yxyy] << "; "
       << "R4[yxyz] = " << R4[yxyz] << endl
       << "R4[yxzx] = " << R4[yxzx] << "; "
       << "R4[yxzy] = " << R4[yxzy] << "; "
       << "R4[yxzz] = " << R4[yxzz] << endl
       << endl;
  cout << "R4[yyxx] = " << R4[yyxx] << "; "
       << "R4[yyxy] = " << R4[yyxy] << "; "
       << "R4[yyxz] = " << R4[yyxz] << endl
       << "R4[yyyx] = " << R4[yyyx] << "; "
       << "R4[yyyy] = " << R4[yyyy] << "; "
       << "R4[yyyz] = " << R4[yyyz] << endl
       << "R4[yyzx] = " << R4[yyzx] << "; "
       << "R4[yyzy] = " << R4[yyzy] << "; "
       << "R4[yyzz] = " << R4[yyzz] << endl
       << endl;
  cout << "R4[yzxx] = " << R4[yzxx] << "; "
       << "R4[yzxy] = " << R4[yzxy] << "; "
       << "R4[yzxz] = " << R4[yzxz] << endl
       << "R4[yzyx] = " << R4[yzyx] << "; "
       << "R4[yzyy] = " << R4[yzyy] << "; "
       << "R4[yzyz] = " << R4[yzyz] << endl
       << "R4[yzzx] = " << R4[yzzx] << "; "
       << "R4[yzzy] = " << R4[yzzy] << "; "
       << "R4[yzzz] = " << R4[yzzz] << endl
       << endl;
  cout << endl;

  cout << "R4[zxxx] = " << R4[zxxx] << "; "
       << "R4[zxxy] = " << R4[zxxy] << "; "
       << "R4[zxxz] = " << R4[zxxz] << endl
       << "R4[zxyx] = " << R4[zxyx] << "; "
       << "R4[zxyy] = " << R4[zxyy] << "; "
       << "R4[zxyz] = " << R4[zxyz] << endl
       << "R4[zxzx] = " << R4[zxzx] << "; "
       << "R4[zxzy] = " << R4[zxzy] << "; "
       << "R4[zxzz] = " << R4[zxzz] << endl
       << endl;
  cout << "R4[zyxx] = " << R4[zyxx] << "; "
       << "R4[zyxy] = " << R4[zyxy] << "; "
       << "R4[zyxz] = " << R4[zyxz] << endl
       << "R4[zyyx] = " << R4[zyyx] << "; "
       << "R4[zyyy] = " << R4[zyyy] << "; "
       << "R4[zyyz] = " << R4[zyyz] << endl
       << "R4[zyzx] = " << R4[zyzx] << "; "
       << "R4[zyzy] = " << R4[zyzy] << "; "
       << "R4[zyzz] = " << R4[zyzz] << endl
       << endl;
  cout << "R4[zzxx] = " << R4[zzxx] << "; "
       << "R4[zzxy] = " << R4[zzxy] << "; "
       << "R4[zzxz] = " << R4[zzxz] << endl
       << "R4[zzyx] = " << R4[zzyx] << "; "
       << "R4[zzyy] = " << R4[zzyy] << "; "
       << "R4[zzyz] = " << R4[zzyz] << endl
       << "R4[zzzx] = " << R4[zzzx] << "; "
       << "R4[zzzy] = " << R4[zzzy] << "; "
       << "R4[zzzz] = " << R4[zzzz] << endl
       << endl;

  cout << "--- Symmetric tensor, rank 4, dimension 2: ---" << endl;
  using sym2_tensor4_t = tensor_u<double, symmetry_type::symmetric, 2, 2, 2, 2>;
  sym2_tensor4_t S4{0};
  S4[xxxx] = 1111;
  S4[yyyy] = 2222;
  S4(0, 0, 0, 1) = 1112;
  S4(0, 0, 1, 1) = 1122;
  S4(1, 0, 1, 1) = 1222;
  cout << "size of S4: " << S4.size() << endl;
  cout << "tensor S4: " << endl << S4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: " << endl;
  cout << "S4[xxxx] = " << S4[xxxx] << "; S4[xxxy] = " << S4[xxxy] << "; "
       << "S4[xxyx] = " << S4[xxyx] << "; S4[xxyy] = " << S4[xxyy] << endl;
  cout << "S4[xyxx] = " << S4[xyxx] << "; S4[xyxy] = " << S4[xyxy] << "; "
       << "S4[xyyx] = " << S4[xyyx] << "; S4[xyyy] = " << S4[xyyy] << endl;
  cout << "S4[yxxx] = " << S4[yxxx] << "; S4[yxxy] = " << S4[yxxy] << "; "
       << "S4[yxyx] = " << S4[yxyx] << "; S4[yxyy] = " << S4[yxyy] << endl;
  cout << "S4[yyxx] = " << S4[yyxx] << "; S4[yyxy] = " << S4[yyxy] << "; "
       << "S4[yyyx] = " << S4[yyyx] << "; S4[yyyy] = " << S4[yyyy] << endl;

  cout << "--- Generic tensor, rank 4, dimension 3: ---" << endl;
  using sym3_tensor4_t = tensor_u<double, symmetry_type::symmetric, 3, 3, 3, 3>;
  sym3_tensor4_t Z4{0};
  Z4[xxxx] = 1111, Z4[yyyy] = 2222;
  Z4[zzzz] = 3333;
  Z4[xxxy] = 1112, Z4[xxyy] = 1122;
  Z4[xyyy] = 1222;
  Z4[xxxz] = 1113, Z4[xxzz] = 1133;
  Z4[xzzz] = 1333;
  Z4[yyyz] = 2223, Z4[yyzz] = 2233;
  Z4[yzzz] = 2333;
  Z4[xyzz] = 1233, Z4[xyyz] = 1223;
  Z4[xxyz] = 1123;

  cout << "size of Z4: " << Z4.size() << endl;
  cout << "tensor Z4: " << endl << Z4 << endl;
  cout << "accessing via the [xyxx.. etc.] operator: " << endl;
  cout << "Z4[xxxx] = " << Z4[xxxx] << "; "
       << "Z4[xxxy] = " << Z4[xxxy] << "; "
       << "Z4[xxxz] = " << Z4[xxxz] << endl
       << "Z4[xxyx] = " << Z4[xxyx] << "; "
       << "Z4[xxyy] = " << Z4[xxyy] << "; "
       << "Z4[xxyz] = " << Z4[xxyz] << endl
       << "Z4[xxzx] = " << Z4[xxzx] << "; "
       << "Z4[xxzy] = " << Z4[xxzy] << "; "
       << "Z4[xxzz] = " << Z4[xxzz] << endl
       << endl;
  cout << "Z4[xyxx] = " << Z4[xyxx] << "; "
       << "Z4[xyxy] = " << Z4[xyxy] << "; "
       << "Z4[xyxz] = " << Z4[xyxz] << endl
       << "Z4[xyyx] = " << Z4[xyyx] << "; "
       << "Z4[xyyy] = " << Z4[xyyy] << "; "
       << "Z4[xyyz] = " << Z4[xyyz] << endl
       << "Z4[xyzx] = " << Z4[xyzx] << "; "
       << "Z4[xyzy] = " << Z4[xyzy] << "; "
       << "Z4[xyzz] = " << Z4[xyzz] << endl
       << endl;
  cout << "Z4[xzxx] = " << Z4[xzxx] << "; "
       << "Z4[xzxy] = " << Z4[xzxy] << "; "
       << "Z4[xzxz] = " << Z4[xzxz] << endl
       << "Z4[xzyx] = " << Z4[xzyx] << "; "
       << "Z4[xzyy] = " << Z4[xzyy] << "; "
       << "Z4[xzyz] = " << Z4[xzyz] << endl
       << "Z4[xzzx] = " << Z4[xzzx] << "; "
       << "Z4[xzzy] = " << Z4[xzzy] << "; "
       << "Z4[xzzz] = " << Z4[xzzz] << endl
       << endl;
  cout << endl;

  cout << "Z4[yxxx] = " << Z4[yxxx] << "; "
       << "Z4[yxxy] = " << Z4[yxxy] << "; "
       << "Z4[yxxz] = " << Z4[yxxz] << endl
       << "Z4[yxyx] = " << Z4[yxyx] << "; "
       << "Z4[yxyy] = " << Z4[yxyy] << "; "
       << "Z4[yxyz] = " << Z4[yxyz] << endl
       << "Z4[yxzx] = " << Z4[yxzx] << "; "
       << "Z4[yxzy] = " << Z4[yxzy] << "; "
       << "Z4[yxzz] = " << Z4[yxzz] << endl
       << endl;
  cout << "Z4[yyxx] = " << Z4[yyxx] << "; "
       << "Z4[yyxy] = " << Z4[yyxy] << "; "
       << "Z4[yyxz] = " << Z4[yyxz] << endl
       << "Z4[yyyx] = " << Z4[yyyx] << "; "
       << "Z4[yyyy] = " << Z4[yyyy] << "; "
       << "Z4[yyyz] = " << Z4[yyyz] << endl
       << "Z4[yyzx] = " << Z4[yyzx] << "; "
       << "Z4[yyzy] = " << Z4[yyzy] << "; "
       << "Z4[yyzz] = " << Z4[yyzz] << endl
       << endl;
  cout << "Z4[yzxx] = " << Z4[yzxx] << "; "
       << "Z4[yzxy] = " << Z4[yzxy] << "; "
       << "Z4[yzxz] = " << Z4[yzxz] << endl
       << "Z4[yzyx] = " << Z4[yzyx] << "; "
       << "Z4[yzyy] = " << Z4[yzyy] << "; "
       << "Z4[yzyz] = " << Z4[yzyz] << endl
       << "Z4[yzzx] = " << Z4[yzzx] << "; "
       << "Z4[yzzy] = " << Z4[yzzy] << "; "
       << "Z4[yzzz] = " << Z4[yzzz] << endl
       << endl;
  cout << endl;

  cout << "Z4[zxxx] = " << Z4[zxxx] << "; "
       << "Z4[zxxy] = " << Z4[zxxy] << "; "
       << "Z4[zxxz] = " << Z4[zxxz] << endl
       << "Z4[zxyx] = " << Z4[zxyx] << "; "
       << "Z4[zxyy] = " << Z4[zxyy] << "; "
       << "Z4[zxyz] = " << Z4[zxyz] << endl
       << "Z4[zxzx] = " << Z4[zxzx] << "; "
       << "Z4[zxzy] = " << Z4[zxzy] << "; "
       << "Z4[zxzz] = " << Z4[zxzz] << endl
       << endl;
  cout << "Z4[zyxx] = " << Z4[zyxx] << "; "
       << "Z4[zyxy] = " << Z4[zyxy] << "; "
       << "Z4[zyxz] = " << Z4[zyxz] << endl
       << "Z4[zyyx] = " << Z4[zyyx] << "; "
       << "Z4[zyyy] = " << Z4[zyyy] << "; "
       << "Z4[zyyz] = " << Z4[zyyz] << endl
       << "Z4[zyzx] = " << Z4[zyzx] << "; "
       << "Z4[zyzy] = " << Z4[zyzy] << "; "
       << "Z4[zyzz] = " << Z4[zyzz] << endl
       << endl;
  cout << "Z4[zzxx] = " << Z4[zzxx] << "; "
       << "Z4[zzxy] = " << Z4[zzxy] << "; "
       << "Z4[zzxz] = " << Z4[zzxz] << endl
       << "Z4[zzyx] = " << Z4[zzyx] << "; "
       << "Z4[zzyy] = " << Z4[zzyy] << "; "
       << "Z4[zzyz] = " << Z4[zzyz] << endl
       << "Z4[zzzx] = " << Z4[zzzx] << "; "
       << "Z4[zzzy] = " << Z4[zzzy] << "; "
       << "Z4[zzzz] = " << Z4[zzzz] << endl
       << endl;

  // cout << "--- Symmetric tensor, rank 3, dimension 4: ---" << endl;
  // using sym4_tensor3_t = tensor_u<double, symmetry_type::symmetric,4,4,4>;
  // for (int i=0;i<4; ++i) {
  //  for (int j=0;j<4; ++j) {
  //    for (int k=0;k<4; ++k) {
  //       cout << sym3_tensor4_t::multiindex(i,j,k) << ", ";
  //    }
  //    cout << endl;
  //  }
  //  cout << endl;
  //}

  // cout << "--- Symmetric tensor, rank 4, dimension 3: ---" << endl;
  // using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,3,3,3,3>;
  // for (int m=0;m<3; ++m) {
  //  for (int i=0;i<3; ++i) {
  //    for (int j=0;j<3; ++j) {
  //      for (int k=0;k<3; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << " -1,"<< endl;
  //    }
  //    cout << "-1, -1, -1, -1," << endl
  //         << endl;
  //  }
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //}

  // cout << "--- Symmetric tensor, rank 4, dimension 2: ---" << endl;
  // using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,2,2,2,2>;
  // for (int m=0;m<2; ++m) {
  //  for (int i=0;i<2; ++i) {
  //    for (int j=0;j<2; ++j) {
  //      for (int k=0;k<2; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << " -1, -1,"<< endl;
  //    }
  //    cout << "-1, -1, -1, -1," << endl
  //         << "-1, -1, -1, -1," << endl
  //         << endl;
  //  }
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //  cout << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << "-1, -1, -1, -1," << endl
  //       << endl;
  //}

  // cout << "--- Symmetric tensor, rank 4, dimension 4: ---" << endl;
  // using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric,4,4,4,4>;
  // for (int m=0;m<4; ++m) {
  //  for (int i=0;i<4; ++i) {
  //    for (int j=0;j<4; ++j) {
  //      for (int k=0;k<4; ++k) {
  //         cout << sym_tensor4_t::multiindex(m,i,j,k) << ", ";
  //      }
  //      cout << endl;
  //    }
  //    cout << endl;
  //  }
  //  cout << endl;
  //}

  // cout << "--- Higher-rank symmetric tensor: ---" << endl;
  // using sym_tensor4_t = tensor_u<double, symmetry_type::symmetric, 3,3,3,3>;
  // sym_tensor4_t Q3{0};
  // cout << "--- Q3 size = " << Q3.size() << endl;
  // Q3(0,0,0,0) = 3333;
  // Q3(0,1,0,1) = 12;
  // Q3(1,2,1,0) = 31;
  // for (int m=0;m<3;++m) {
  //  for (int i=0;i<3; ++i) {
  //    for (int j=0;j<3; ++j) {
  //      for (int k=0;k<3; ++k) {
  //         cout << Q3(m,i,j,k) << " ";
  //      }
  //      cout << endl;
  //    }
  //    cout << " ----------- " << endl;
  //  }
  //  cout << "=============" << endl;
  //}
}