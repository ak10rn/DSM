#include <bits/stdc++.h>
#define pb push_back
#define vd vector<double>
#define vvd vector<vd>
#define vi vector<int>
#define vvi vector<vi>
#define vdd vector<pair<double, double>>
#define vvT vector<vector<T>>
using namespace std;


auto display(const vvd& A, const vi& B = {}, bool top = true) {
  cout << '\n';
  bool flag = false;
  string tmp = "---------------", s = "--------";
  if (B.size() and top) {
    for (auto cell : B) {
      cout << setprecision(10) << setw(15) << setfill(' ') << cell << ' ';
      s += tmp;
    }
    cout << "\n        " << s;
  }
  if (B.size()) flag = true;
  cout << '\n';
  auto it = B.begin();
  for (auto& row : A) {
    for (auto cell : row)
      cout << setprecision(10) << setw(15) << setfill(' ') << cell << ' ';
    if (flag) {
      cout << setprecision(10) << setw(15) << setfill(' ') << "|  " << *it
           << "      ";
      it++;
    }
    cout << '\n';
  }
  cout << '\n';
}

auto input() {
  cout << "Enter matrix and press extra enter when done!!\n";
  string s;
  getline(cin, s);
  vvd ans;
  while (s != "") {
    vd row;
    stringstream ss(s);
    double cell;
    while (ss >> cell) row.pb(cell);
    ans.pb(row);
    getline(cin, s);
  }
  return ans;
}

auto getCoFactor(const vvd& A, const int& r, const int& c) {
  int n = A.size();
  vvd ans;
  for (int i = 0; i < n; i++) {
    if (i == r) continue;
    vd v;
    for (int j = 0; j < n; j++) {
      if (j == c) continue;
      v.pb(A[i][j]);
    }
    ans.pb(v);
  }
  return ans;
}

auto determinant(const vvd& A) {
  int n = A.size();
  if (n == 1) return A[0][0];
  double ans = 0;
  for (int i = 0; i < n; i++) {
    vvd co = getCoFactor(A, 0, i);
    ans += pow(-1, i) * A[0][i] * determinant(co);
  }
  return ans;
}

auto inverseMatrix(const vvd& A) {
  if (A.size() == 1) {
    vvd ans = {{1.0 / A[0][0]}};
    return ans;
  }
  double det = determinant(A);
  if (det == 0) {
    cout << "Non-Invertible Matrix";
    exit(1);
  }
  int n = A.size();
  vvd ans(n, vd(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      ans[j][i] = pow(-1, i + j) * determinant(getCoFactor(A, i, j)) / det;
    }
  }
  return ans;
}

auto multiplyMatrix(const vvd& A, const vvd& B) {
  int r1 = A.size(), c1 = A[0].size();
  int r2 = B.size(), c2 = B[0].size();
  if (r2 != c1) {
    cout << "Multiplication Not possible!!!";
    exit(1);
  }
  vvd ans(r1, vd(c2, 0.0));
  for (int i = 0; i < r1; i++) {
    for (int j = 0; j < c2; j++) {
      for (int k = 0; k < c1; k++) {
        ans[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return ans;
}

auto multiplyConstant(vvd& A, const auto& mul) {
  for (auto& i : A)
    for (auto& j : i) j *= mul;
}

auto divideConstant(vvd& A, const auto& den) {
  for (auto& i : A)
    for (auto& j : i) j /= den;
}

auto add(const vvd& A, const vvd& B, bool isAdding = true) {
  int n = A.size(), m = A[0].size();
  vvd ans(n, vd(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (isAdding)
        ans[i][j] = A[i][j] + B[i][j];
      else
        ans[i][j] = A[i][j] - B[i][j];
    }
  }
  return ans;
}

auto transpose(const vvd& A) {
  vvd ans(A[0].size(), vd(A.size()));
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < A[0].size(); j++) {
      ans[j][i] = A[i][j];
    }
  }
  return ans;
}

auto gk(const vector<vvd>& allEKG, const vvi& designation) {
  int n = 0;
  for (const auto& i : designation) n = max(n, *max_element(begin(i), end(i)));
  vvd GK(n, vd(n, 0));
  for (int i = 0; i < allEKG.size(); i++) {
    vi des = designation[i];
    vvd ekg = allEKG[i];
    for (int j = 0; j < ekg.size(); j++) {
      for (int k = 0; k < ekg[0].size(); k++) {
        GK[des[j] - 1][des[k] - 1] += ekg[j][k];
      }
    }
  }
  return GK;
}

auto getEKG(const vvd& ek, const vvd& T, const auto& l = 1) {
  auto ans = multiplyMatrix(transpose(T), multiplyMatrix(ek, T));
  if (l != 1) divideConstant(ans, l);
  return ans;
}

auto openOutputFile(const string& filename) {
  string path = filesystem::current_path().string(), command = "cd ";
  for (char c : path) command += c, c == '\\' ? command += c : "";
  command += " && " + filename;
  system(command.c_str());
}

auto Truss() {
  cout << "\nBe Careful about units :) \n\n";
  cout << "Enter number of nodes: ";
  int nodes, members;
  cin >> nodes;
  cout << "Enter number of members: ";
  cin >> members;
  vvi designation;
  vector<vvd> allEKG, allTransformation;
  vd memberWiseLength, memberWiseAxialRigidity;
  vvd memberWiseLambda;
  vvd qo;
  vd memberwiseErrorTermdl_l;

  for (int i = 0; i < members; i++) {
    cout << "Enter Nx Ny Fx Fy for member " << i + 1 << " : ";
    vi d(4);
    for (auto& i : d) cin >> i;
    designation.pb(d);
    cout << "Is any support is inclined (y/n): ";
    string s;
    cin >> s;
    double x, y, x1, y1;
    if (s == "Y" or s == "y") {
      cout << "Enter lambdax lambday lambdax'' lambday'' for member " << i + 1
           << " : ";
      cin >> x >> y >> x1 >> y1;
    } else {
      cout << "Enter lambdax lambday for member " << i + 1 << " : ";
      cin >> x >> y;
      x1 = x, y1 = y;
    }
    memberWiseLambda.pb({x, y, x1, y1});
    vvd T = {{x, y, 0, 0}, {0, 0, x1, y1}};
    allTransformation.pb(T);
    cout << "Enter length of member " << i + 1 << " : ";
    double l;
    cin >> l;
    memberWiseLength.pb(l);
    cout << "Enter Axial Rigidity of member " << i + 1 << " : ";
    double ae;
    cin >> ae;
    memberWiseAxialRigidity.pb(ae);
    cout << "Is it contains Thermal Loading or Fabrivation Error (y/n): ";
    cin >> s;
    double dl_l = 0;
    if (s == "Y" or s == "y") {
      cout << "Enter value of (delta L)/L or  (alpha delta T): ";
      cin >> dl_l;
    }
    memberwiseErrorTermdl_l.pb(dl_l);
    qo.pb({ae * dl_l * x, ae * dl_l * y, -ae * dl_l * x1, -ae * dl_l * y1});
    vvd kg = {{1, -1}, {-1, 1}};
    multiplyConstant(kg, ae);
    allEKG.pb(getEKG(kg, T, l));
    cout << '\n';
  }

  vvd Qo(2 * nodes, vd(1, 0));
  for (int i = 0; i < members; i++) {
    auto d = designation[i];
    for (int j = 0; j < d.size(); j++) {
      Qo[d[j] - 1][0] += qo[i][j];
    }
  }
  vvd GK = gk(allEKG, designation);
  vvd D(2 * nodes), Q(2 * nodes), Q_known;
  vector<bool> D_isKnown(2 * nodes, true), Q_isKnown(2 * nodes, true);
  cout << "Enter known forces sequentially, if not known enter '.' : ";
  for (int i = 0; i < 2 * nodes; i++) {
    string s;
    cin >> s;
    if (s == ".")
      Q_isKnown[i] = false;
    else
      Q[i] = {stod(s)};
  }
  cout << "Enter known displacement sequentially, if not known enter '.' : ";
  for (int i = 0; i < 2 * nodes; i++) {
    string s;
    cin >> s;
    if (s == ".")
      D_isKnown[i] = false;
    else
      D[i] = {stod(s)};
  }
  cout << "\n\n";
  vvd K11, K12, K21, K22, Qko;
  for (int i = 0; i < 2 * nodes; i++) {
    if (!D_isKnown[i]) {
      vd temp1, temp2;
      Qko.pb({Qo[i]});
      for (int j = 0; j < 2 * nodes; j++) {
        if (!D_isKnown[j])
          temp1.pb(GK[i][j]);
        else
          temp2.pb(GK[i][j]);
      }
      if (temp1.size()) K11.pb(temp1);
      if (temp2.size()) K12.pb(temp2);
    } else {
      vd temp1, temp2;
      for (int j = 0; j < 2 * nodes; j++) {
        if (!D_isKnown[j])
          temp1.pb(GK[i][j]);
        else
          temp2.pb(GK[i][j]);
      }
      if (temp1.size()) K21.pb(temp1);
      if (temp2.size()) K22.pb(temp2);
    }
  }
  for (int i = 0; i < 2 * nodes; i++) {
    if (Q_isKnown[i]) Q_known.pb({Q[i]});
  }
  vvd D_Known;
  for (int i = 0; i < 2 * nodes; i++) {
    if (D_isKnown[i]) D_Known.pb(D[i]);
  }
  auto K12xD_Known = multiplyMatrix(K12, D_Known);
  auto kachra = add(K12xD_Known, Qko);
  auto D_unknown =
      multiplyMatrix(inverseMatrix(K11), add(Q_known, kachra, false));
  int pos = 0;
  for (int i = 0; i < 2 * nodes; i++) {
    if (!D_isKnown[i]) D[i] = D_unknown[pos++], D_isKnown[i] = true;
  }
  Q = multiplyMatrix(GK, D);
  for (int i = 0; i < 2 * nodes; i++) Q_isKnown[i] = true;

  vd localMemberForces;
  for (int i = 0; i < members; i++) {
    double x = memberWiseLambda[i][0];
    double y = memberWiseLambda[i][1];
    double x1 = memberWiseLambda[i][2];
    double y1 = memberWiseLambda[i][3];
    double l = memberWiseLength[i];
    double ae = memberWiseAxialRigidity[i];
    vvd temp;
    for (int j : designation[i]) temp.pb(D[j - 1]);
    temp = multiplyMatrix({{-x, -y, x1, y1}}, temp);
    multiplyConstant(temp, 1.0 * ae / l);
    localMemberForces.pb(temp[0][0] - ae * memberwiseErrorTermdl_l[i]);
  }
  vi allDOF(2 * nodes);
  iota(allDOF.begin(), allDOF.end(), 1);

  string fname = "DSM_Output.txt";
  auto fp = freopen(fname.c_str(), "w", stdout);
  cout << "Memberwise EKG\n";
  for (int i = 0; i < members; i++) {
    cout << "Member " << i + 1 << '\n';
    display(allEKG[i], designation[i]);
  }
  cout << "Global Stiffness Matrix\n";
  display(GK, allDOF);
  cout << "Global Force Matrix\n";
  display(Q, allDOF, false);
  cout << "Global Displacement Matrix\n";
  display(D, allDOF, false);
  cout << "Memberwise Local Displacement Matrix\n";
  for (int i = 0; i < members; i++) {
    cout << "Member " << i + 1 << '\n';
    vvd D1;
    for (int j : designation[i]) D1.pb(D[j - 1]);
    display(multiplyMatrix(allTransformation[i], D1));
  }
  cout << "Local Member Forces (-ve for tensile and +ve for compressive)\n";
  display(transpose({localMemberForces}));

  fclose(fp);
  openOutputFile(fname);
}

auto Beam() {
  cout << "\nBe Careful about units :) \n\n";
  cout << "Enter number of nodes: ";
  int nodes, members;
  cin >> nodes;
  cout << "Enter number of members: ";
  cin >> members;
  vvi designation;
  vector<vvd> allEKG;
  vd memberWiseLength, memberWiseFlexuralRigidity;
  vvd qo;
  for (int i = 0; i < members; i++) {
    cout << "Enter Ny Nz Fy Fz for member " << i + 1 << " : ";
    vi d(4);
    for (auto& i : d) cin >> i;
    designation.pb(d);
    cout << "Enter length of member " << i + 1 << " : ";
    double l;
    cin >> l;
    memberWiseLength.pb(l);
    cout << "Enter Flexural Rigidity of member " << i + 1 << " : ";
    double ei;
    cin >> ei;
    memberWiseFlexuralRigidity.pb(ei);
    vvd EKG = {{12.0, 6.0 * l, -12.0, 6.0 * l},
               {6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l},
               {-12.0, -6.0 * l, 12.0, -6.0 * l},
               {6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l}};
    multiplyConstant(EKG, 12 * ei / (l * l * l));
    allEKG.pb(EKG);
    cout << "Is this member containes non nodal loads (y/n): ";
    string s;
    cin >> s;
    vd fem(4, 0);
    if (s == "y" or s == "Y") {
      cout << "Enter Fixed End Force Matrix respective to Designation given: ";
      for (double& i0 : fem) cin >> i0;
    }
    qo.pb(fem);
    cout << '\n';
  }
  vvd GK = gk(allEKG, designation);

  vvd D(2 * nodes, vd(1, 0)), Q(2 * nodes, vd(1, 0)), Q_known;
  vector<bool> D_isKnown(2 * nodes, true), Q_isKnown(2 * nodes, true);
  cout << "Enter known forces sequentially (do not include FEM), if not known "
          "enter '.' : ";
  for (int i = 0; i < 2 * nodes; i++) {
    string s;
    cin >> s;
    if (s == ".")
      Q_isKnown[i] = false;
    else
      Q[i] = {stod(s)};
  }
  cout << "Enter known displacement sequentially, if not known enter '.' : ";
  for (int i = 0; i < 2 * nodes; i++) {
    string s;
    cin >> s;
    if (s == ".")
      D_isKnown[i] = false;
    else
      D[i] = {stod(s)};
  }
  for (int i = 0; i < members; i++) {
    int pos = 0;
    for (auto j : designation[i]) {
      if (Q_isKnown[j - 1]) Q[j - 1][0] -= qo[i][pos];
      pos++;
    }
  }  // 2 3 2 5 6 2 1 6 1 y 27 27 27 -27 2 1 3 4 4 1 y 12 8 12 -8 0 . . . . . . 0 0 0 0 0
  // display(Q);
  cout << "\n\n";
  vvd K11, K12, K21, K22, Qko;
  for (int i = 0; i < 2 * nodes; i++) {
    if (!D_isKnown[i]) {
      vd temp1, temp2;
      for (int j = 0; j < 2 * nodes; j++) {
        if (!D_isKnown[j])
          temp1.pb(GK[i][j]);
        else
          temp2.pb(GK[i][j]);
      }
      if (temp1.size()) K11.pb(temp1);
      if (temp2.size()) K12.pb(temp2);
    } else {
      vd temp1, temp2;
      for (int j = 0; j < 2 * nodes; j++) {
        if (!D_isKnown[j])
          temp1.pb(GK[i][j]);
        else
          temp2.pb(GK[i][j]);
      }
      if (temp1.size()) K21.pb(temp1);
      if (temp2.size()) K22.pb(temp2);
    }
  }
  for (int i = 0; i < 2 * nodes; i++) {
    if (Q_isKnown[i]) Q_known.pb({Q[i]});
  }
  vvd D_Known;
  for (int i = 0; i < 2 * nodes; i++) {
    if (D_isKnown[i]) D_Known.pb(D[i]);
  }
  auto K12xD_Known = multiplyMatrix(K12, D_Known);
  auto tempppp = add(Q_known, K12xD_Known, false);  //
  auto D_unknown = multiplyMatrix(inverseMatrix(K11), tempppp);
  // display(K11);
  // display(K12);
  // display(K21);
  // display(K22);
  // display(D_unknown);
  int pos = 0;
  for (int i = 0; i < 2 * nodes; i++) {
    if (!D_isKnown[i]) D[i] = D_unknown[pos++], D_isKnown[i] = true;
  }
  Q = multiplyMatrix(GK, D);
  for (int i = 0; i < 2 * nodes; i++) Q_isKnown[i] = true;

  vector<vvd> memberWiseReactions;
  for (auto& i : Q) i = {0};
  for (int i = 0; i < members; i++) {
    vvd Disp;
    for (int j : designation[i]) Disp.pb(D[j - 1]);
    auto temp = add(multiplyMatrix(allEKG[i], Disp), transpose({qo[i]}));
    memberWiseReactions.pb(temp);
    int pos11 = 0;
    for (int j : designation[i]) Q[j - 1][0] += temp[pos11++][0];
  }

  vi allDOF(2 * nodes);
  iota(allDOF.begin(), allDOF.end(), 1);

  string fname = "DSM_Output.txt";
  auto fp = freopen(fname.c_str(), "w", stdout);
  cout << "Memberwise EKG\n";
  for (int i = 0; i < members; i++) {
    cout << "Member " << i + 1 << '\n';
    display(allEKG[i], designation[i]);
  }
  cout << "Global Stiffness Matrix\n";
  display(GK, allDOF);
  cout << "Global Force Matrix\n";
  display(Q, allDOF, false);
  cout << "Global Displacement Matrix\n";
  display(D, allDOF, false);
  cout << "Memberwise Reactions\n";
  for (int i = 0; i < members; i++) {
    cout << "Member " << i + 1 << '\n';
    display(memberWiseReactions[i], designation[i], false);
  }

  fclose(fp);
  openOutputFile(fname);
}

auto Frame() {
  char normal[] = {0x1b, '[', '0', ';', '3', '9', 'm', 0},
       Bred[] = {0x1b, '[', '1', ';', '3', '1', 'm', 0},
       Emoji[] = {0x1b, '[', '1', ';', '3', '1', 'm', 1};
  string s = "\n\nBacche ki Jaan Loge Kya  ";
  cout << Bred << s << Emoji << normal << "\n\n\n\n";
  // system("shutdown -h now");

  // cout << "\nBe Careful about units :) \n\n";
  // cout << "Enter number of nodes: ";
  // int nodes, members;
  // cin >> nodes;
  // cout << "Enter number of members: ";
  // cin >> members;
  // vvi designation;
  // vector<vvd> allEKG;
  // vd memberWiseLength;
  // vdd memberWiseLambda, memberWiseAxialRigidity;
  // for (int i = 0; i < members; i++) {
  //   cout << "Enter Nx Ny Nz Fx Fy Fz for member " << i + 1 << " : ";
  //   vi d(6);
  //   for (auto& i : d) cin >> i;
  //   designation.pb(d);
  //   cout << "Enter lambdax lambday for member " << i + 1 << " : ";
  //   double x, y;
  //   cin >> x >> y;
  //   memberWiseLambda.pb({x, y});
  //   cout << "Enter length of member " << i + 1 << " : ";
  //   double l;
  //   cin >> l;
  //   memberWiseLength.pb(l);
  //   cout << "Enter AE EI of member " << i + 1 << " : ";
  //   double ae, ei;
  //   cin >> ae >> ei;
  //   memberWiseAxialRigidity.pb({ae, ei});
  //   vvd T = {{x, -y, 0, 0, 0, 0}, {y, x, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0},
  //            {0, 0, 0, x, -y, 0}, {0, 0, 0, y, x, 0}, {0, 0, 0, 0, 0, 1}};
  //   vvd kg = {{ae / l, 0, 0, -ae / l, 0, 0},
  //             {0, 12.0 * ei / (l * l * l), 6.0 * ei / (l * l), 0,
  //              -12.0 * ei / (l * l * l), 6.0 * ei / (l * l)},
  //             {0, 6.0 * ei / (l * l), 4.0 * ei / l, 0, -6.0 * ei / (l * l),
  //              2.0 * ei / l},
  //             {-ae / l, 0, 0, ae / l, 0, 0},
  //             {0, -12.0 * ei / (l * l * l), -6.0 * ei / (l * l), 0,
  //              12.0 * ei / (l * l * l), -6.0 * ei / (l * l)},
  //             {0, 6.0 * ei / (l * l), 2.0 * ei / l, 0, -6.0 * ei / (l * l),
  //              4.0 * ei / l}};
  //   allEKG.pb(getEKG(kg, T, 1));
  //   cout << '\n';
  // }
  // vvd GK = gk(allEKG, designation);
  // vvd D(3 * nodes), Q(3 * nodes), Q_known, K_partitioned;
  // vector<bool> D_isKnown(3 * nodes, true), Q_isKnown(3 * nodes, true);
  // cout << "Enter known forces sequentially, if not known enter '.' : ";
  // for (int i = 0; i < 3 * nodes; i++) {
  //   string s;
  //   cin >> s;
  //   if (s == ".")
  //     Q_isKnown[i] = false;
  //   else
  //     Q[i] = {stod(s)};
  // }
  // cout << "Enter known displacement sequentially, if not known enter '.' : ";
  // for (int i = 0; i < 3 * nodes; i++) {
  //   string s;
  //   cin >> s;
  //   if (s == ".")
  //     D_isKnown[i] = false;
  //   else
  //     D[i] = {stod(s)};
  // }
  // cout << "\n\n";
  // for (int i = 0; i < 3 * nodes; i++) {
  //   if (D_isKnown[i]) continue;
  //   vd temp;
  //   for (int j = 0; j < 3 * nodes; j++) {
  //     if (!D_isKnown[j]) temp.pb(GK[i][j]);
  //   }
  //   K_partitioned.pb(temp);
  // }
  // for (int i = 0; i < 3 * nodes; i++) {
  //   if (Q_isKnown[i]) Q_known.pb({Q[i]});
  // }
  // auto D_unknown = multiplyMatrix(inverseMatrix(K_partitioned), Q_known);
  // int pos = 0;
  // for (int i = 0; i < 3 * nodes; i++) {
  //   if (!D_isKnown[i]) D[i] = D_unknown[pos++], D_isKnown[i] = true;
  // }
  // Q = multiplyMatrix(GK, D);
  // for (int i = 0; i < 3 * nodes; i++) Q_isKnown[i] = true;

  // // check

  // // vd localMemberForces;

  // // for (int i = 0; i < members; i++) {
  // //   auto [x, y] = memberWiseLambda[i];
  // //   double l = memberWiseLength[i];
  // //   double ae = memberWiseAxialRigidity[i];
  // //   vvd temp;
  // //   for (int j : designation[i]) temp.pb(D[j - 1]);
  // //   temp = multiplyMatrix({{x, -y, 0, 0, 0, 0},
  // //                          {y, x, 0, 0, 0, 0},
  // //                          {0, 0, 1, 0, 0, 0},
  // //                          {0, 0, 0, x, -y, 0},
  // //                          {0, 0, 0, y, x, 0},
  // //                          {0, 0, 0, 0, 0, 1}},
  // //                         temp);
  // //   multiplyConstant(temp, 1.0 * ae / l);
  // //   localMemberForces.pb(temp[0][0]);
  // // }

  // //

  // string fname = "DSM_Output.txt";
  // auto fp = freopen(fname.c_str(), "w", stdout);
  // cout << "Global Stiffness Matrix\n";
  // display(GK);
  // cout << "Global Force Matrix\n";
  // display(Q);
  // cout << "Global Displacement Matrix\n";
  // display(D);
  // cout << "Memberwise EKG\n";
  // for (int i = 0; i < members; i++) {
  //   cout << "Member " << i + 1 << '\n';
  //   display(allEKG[i]);
  // }

  // fclose(fp);
  // openOutputFile(fname);
}

void getEKGOnly() {
  cout << "1 Truss\n2 Beam\n3 Frame\nEnter your choice: ";
  int choice;
  cin >> choice;
  system("CLS");
  switch (choice) {
    case 1: {
      double x, y, l, x1, y1, ae;
      cout << "Enter labmdax lambday lambdax'' lambday'' length: ";
      cin >> x >> y >> x1 >> y1 >> l;
      cout << "Enter Axial Rigidity: ";
      cin >> ae;
      vvd T = {{x, y, 0, 0}, {0, 0, x1, y1}};
      vvd kg = {{ae, -ae}, {-ae, ae}};
      string fname = "DSM_Output.txt";
      auto fp = freopen(fname.c_str(), "w", stdout);
      display(getEKG(kg, T, l));
      fclose(fp);
      openOutputFile(fname);
      break;
    }
    case 2: {
      double l, ei;
      cout << "Enter length and Flexural Rigidity: ";
      cin >> l >> ei;
      vvd EKG = {{12.0, 6.0 * l, -12.0, 6.0 * l},
                 {6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l},
                 {-12.0, -6.0 * l, 12.0, -6.0 * l},
                 {6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l}};
      multiplyConstant(EKG, 12 * ei / (l * l * l));
      display(EKG);
      break;
    }
    case 3: {
      double x, y, l, x1, y1, ae, ei;
      cout << "Enter labmdax lambday lambdax'' lambday'' length: ";
      cin >> x >> y >> x1 >> y1 >> l;
      cout << "Enter Axial Rigidity and Flexural Rigidity: ";
      cin >> ae >> ei;
      vvd T = {{x, -y, 0, 0, 0, 0}, {y, x, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, x, -y, 0}, {0, 0, 0, y, x, 0}, {0, 0, 0, 0, 0, 1}};
      vvd kg = {{ae / l, 0, 0, -ae / l, 0, 0},
                {0, 12.0 * ei / (l * l * l), 6.0 * ei / (l * l), 0,
                 -12.0 * ei / (l * l * l), 6.0 * ei / (l * l)},
                {0, 6.0 * ei / (l * l), 4.0 * ei / l, 0, -6.0 * ei / (l * l),
                 2.0 * ei / l},
                {-ae / l, 0, 0, ae / l, 0, 0},
                {0, -12.0 * ei / (l * l * l), -6.0 * ei / (l * l), 0,
                 12.0 * ei / (l * l * l), -6.0 * ei / (l * l)},
                {0, 6.0 * ei / (l * l), 2.0 * ei / l, 0, -6.0 * ei / (l * l),
                 4.0 * ei / l}};
      display(getEKG(kg, T, 1));
      break;
    }
    default:
      getEKGOnly();
  }
}

void solve() {
  cout << "1 Truss\n2 Beam\n3 Frame\n4 EKG Only\nEnter your choice: ";
  int choice;
  cin >> choice;
  system("CLS");
  switch (choice) {
    case 1:
      Truss();
      break;
    case 2:
      Beam();
      break;
    case 3:
      Frame();
      break;
    case 4:
      getEKGOnly();
      break;
    default:
      solve();
  }
}

auto main() -> int {
  system("CLS");
  solve();
}

// F

// Beam 2 3 2 4 3 5 2 288 15790000 y 24 1152 24 -1152 5 2 6 1 96 14790000 y 6
// 144 6 -144 0 0 . . . . . . 0 0 0 0
