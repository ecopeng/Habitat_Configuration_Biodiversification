/*

Functions

*/

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cmath>
#include "igraph.h"
#include "Appendix_4_poiss_rng.c" // the Poisson random number generator published by Ashby, B., Shaw, A. K., & Kokko, H. (2020, Oikos): https://doi.org/10.1111/oik.06704

int SAMPLE(float p) {
  std::vector<float> Prob = {1 - p, p};
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> dis_distr(Prob.begin(), Prob.end());
  return dis_distr(gen);
}

int SAMPLE(igraph_vector_t* V) {
  std::vector<float> Prob(igraph_vector_size(V));
  std::fill(Prob.begin(), Prob.end(), 1);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> dis_distr(Prob.begin(), Prob.end());
  return VECTOR(*V)[dis_distr(gen)];
}

float Normal(float mean, float sd) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<float> normal{mean, sd};
  return normal(gen);;
}

igraph_matrix_t SUB_MAT(float g1, float g2, int nt) {
  igraph_matrix_t compe;
  igraph_matrix_init(&compe, nt - 1, nt - 1);
  for(int i = 0; i < nt - 1; ++i) {
    for(int j = 0; j < nt - 1; ++j) {
      igraph_matrix_set(&compe, i, j, exp(-g1 * i * i) * exp(-g2 * j * j));
    }
  }
  return compe;
}

igraph_matrix_t UNS_MAT(float g1, float g2, int nt, float w1, float w2) {
  igraph_matrix_t compe;
  igraph_matrix_init(&compe, nt - 1, nt - 1);
  for(int i = 0; i < nt - 1; ++i) {
    for(int j = 0; j < nt - 1; ++j) {
      igraph_matrix_set(&compe, i, j, w1 * exp(-g1 * i * i) + w2 * exp(-g2 * j * j));
    }
  }
  return compe;
}

void GROWTH(igraph_matrix_t* compe, igraph_matrix_t* N, igraph_matrix_t* K, igraph_vector_t* T1, igraph_vector_t* T2, float r, int nt) {
  igraph_vector_t temp;
  igraph_vector_init(&temp, igraph_matrix_ncol(N));
  for(int x = 0; x < igraph_matrix_nrow(N); ++x) {
    igraph_matrix_get_row(N, &temp, x);
    if(igraph_vector_sum(&temp) != 0) {
      float comp;
      for(int i = 0; i < igraph_vector_size(&temp); ++i) {
        if(igraph_vector_e(&temp, i) != 0) {
          comp = 0;
          for(int j = 0; j < igraph_vector_size(&temp); ++j) {
            comp += igraph_matrix_e(compe, abs(igraph_vector_e(T1, i) - igraph_vector_e(T1, j)) == nt - 1 ? 1 : abs(igraph_vector_e(T1, i) - igraph_vector_e(T1, j)), abs(igraph_vector_e(T2, i) - igraph_vector_e(T2, j)) == nt - 1 ? 1 : abs(igraph_vector_e(T2, i) - igraph_vector_e(T2, j))) * igraph_vector_e(&temp, j);
          }
          igraph_matrix_set(N, x, i, poiss_large(std::max<float>(0, 1 + r * (1 - comp / igraph_matrix_e(K, x, i))) * igraph_vector_e(&temp, i)));
        }
      }
    }
  }
  igraph_vector_destroy(&temp);
}

void SPECIATION(igraph_matrix_t* N, igraph_vector_t* T1, igraph_vector_t* T2, igraph_vector_t* T3, float mu, int ck) {
  int n1, n2, n3, n4;
  igraph_vector_t temp;
  igraph_vector_init(&temp, igraph_matrix_ncol(N));
  for(int x = 0; x < igraph_matrix_nrow(N); ++x) {
    igraph_matrix_get_row(N, &temp, x);
    if(igraph_vector_sum(&temp) != 0) {
      for(int i = 0; i < igraph_vector_size(&temp); ++i) {
        if(igraph_vector_e(&temp, i) != 0) {
          int t_1_L = igraph_vector_e(T1, i) - 1 >= 0 ? igraph_vector_e(T1, i) - 1 : ck - 1;
          int t_1_R = igraph_vector_e(T1, i) + 1 > ck - 1 ? 0 : igraph_vector_e(T1, i) + 1;
          int t_2_D = igraph_vector_e(T2, i) - 1 >= 0 ? igraph_vector_e(T2, i) - 1 : ck - 1;
          int t_2_U = igraph_vector_e(T2, i) + 1 > ck - 1 ? 0 : igraph_vector_e(T2, i) + 1;
          for(int j = 0; j < igraph_vector_size(&temp); ++j) {
            if(j != i) {
              n1 = poiss_large(mu * igraph_vector_e(&temp, i));
              n2 = poiss_large(mu * igraph_vector_e(&temp, i));
              n3 = poiss_large(mu * igraph_vector_e(&temp, i));
              n4 = poiss_large(mu * igraph_vector_e(&temp, i));
              while(n1 + n2 + n3 + n4 > igraph_vector_e(&temp, i)) {
                n1 = poiss_large(mu * igraph_vector_e(&temp, i));
                n2 = poiss_large(mu * igraph_vector_e(&temp, i));
                n3 = poiss_large(mu * igraph_vector_e(&temp, i));
                n4 = poiss_large(mu * igraph_vector_e(&temp, i));
              }
              if(n1 != 0 && igraph_vector_e(T1, j) == t_1_L && igraph_vector_e(T2, j) == igraph_vector_e(T2, i) && igraph_vector_e(T3, j) == x) {
                igraph_matrix_set(N, x, i, igraph_vector_e(&temp, i) - n1);
                igraph_matrix_set(N, x, j, igraph_vector_e(&temp, j) + n1);
              }
              if(n2 != 0 && igraph_vector_e(T1, j) == t_1_R && igraph_vector_e(T2, j) == igraph_vector_e(T2, i) && igraph_vector_e(T3, j) == x) {
                igraph_matrix_set(N, x, i, igraph_vector_e(&temp, i) - n2);
                igraph_matrix_set(N, x, j, igraph_vector_e(&temp, j) + n2);
              }
              if(n3 != 0 && igraph_vector_e(T2, j) == t_2_D && igraph_vector_e(T1, j) == igraph_vector_e(T1, i) && igraph_vector_e(T3, j) == x) {
                igraph_matrix_set(N, x, i, igraph_vector_e(&temp, i) - n3);
                igraph_matrix_set(N, x, j, igraph_vector_e(&temp, j) + n3);
              }
              if(n4 != 0 && igraph_vector_e(T2, j) == t_2_U && igraph_vector_e(T1, j) == igraph_vector_e(T1, i) && igraph_vector_e(T3, j) == x) {
                igraph_matrix_set(N, x, i, igraph_vector_e(&temp, i) - n4);
                igraph_matrix_set(N, x, j, igraph_vector_e(&temp, j) + n4);
              }
            }
          }
        }
      }
    }
  }
  igraph_vector_destroy(&temp);
}

void DISPERSAL(igraph_matrix_t* N, float rho, igraph_matrix_t* habitat) {
  igraph_matrix_t Nc;
  igraph_matrix_copy(&Nc, N);
  igraph_vector_t temp, nD;
  igraph_vector_init(&temp, igraph_matrix_ncol(habitat));
  igraph_vector_init(&nD, igraph_matrix_ncol(habitat));
  for(int i = 0; i < igraph_matrix_nrow(habitat); ++i) {
    igraph_matrix_get_row(habitat, &temp, i);
    igraph_vector_scale(&temp, rho / igraph_vector_sum(&temp));
    for(int k = 0; k < igraph_matrix_ncol(&Nc); ++k) {
      if(igraph_matrix_e(&Nc, i, k) != 0) {
        for(int l = 0; l < igraph_vector_size(&nD); ++l) {
          igraph_vector_set(&nD, l, poiss_large(VECTOR(temp)[l] * igraph_matrix_e(&Nc, i, k)));
        }
        while(igraph_vector_sum(&nD) > igraph_matrix_e(&Nc, i, k)) {
          for(int l = 0; l < igraph_vector_size(&nD); ++l) {
            igraph_vector_set(&nD, l, poiss_large(VECTOR(temp)[l] * igraph_matrix_e(&Nc, i, k)));
          }
        }
        for(int l = 0; l < igraph_vector_size(&nD); ++l) {
          if(VECTOR(nD)[l] != 0) {
            igraph_matrix_set(N, i, k, igraph_matrix_e(N, i, k) - VECTOR(nD)[l]);
            igraph_matrix_set(N, l, k, igraph_matrix_e(N, l, k) + VECTOR(nD)[l]);
          }
        }
      }
    }
  }
  igraph_matrix_destroy(&Nc);
  igraph_vector_destroy(&temp);
  igraph_vector_destroy(&nD);
}

void EXTINCTION(igraph_matrix_t* N, float er) {
  for(int i = 0; i < igraph_matrix_nrow(N); ++i) {
    for(int j = 0; j < igraph_matrix_ncol(N); ++j) {
      if(igraph_matrix_e(N, i, j) != 0 && SAMPLE(er)) {
        igraph_matrix_set(N, i, j, 0);
      }
    }
  }
}

int RICHNESS(igraph_matrix_t* N, int delta) {
  igraph_vector_t res;
  igraph_vector_init(&res, igraph_matrix_ncol(N));
  igraph_matrix_colsum(N, &res);
  int r = 0;
  for(int i = 0; i < igraph_vector_size(&res); ++i) {
    if(VECTOR(res)[i] > delta) {
      r += 1;
    }
  }
  igraph_vector_destroy(&res);
  return r;
}

void Initializaiton(igraph_vector_t* T1, igraph_vector_t* T2, igraph_vector_t* T3, igraph_matrix_t* K, igraph_matrix_t* N, int nt, int np, float mean, float sd) {
  for(int t3 = 0; t3 < np; ++t3) {
    for(int t2 = 0; t2 < nt; ++t2) {
      for(int t1 = 0; t1 < nt; ++t1) {
        igraph_vector_push_back(T3, t3);
        igraph_vector_push_back(T2, t2);
        igraph_vector_push_back(T1, t1);
      }
    }
  }
  igraph_vector_t vec;
  igraph_vector_init_seq(&vec, 0, np - 1);
  int r = SAMPLE(&vec);
  igraph_matrix_set(N, r, (nt * nt - 1) / 2 + r * nt * nt, std::max<float>(Normal(mean, sd), 0));
  igraph_vector_destroy(&vec);
  for(int i = 0; i < igraph_matrix_nrow(K); ++i) {
    for(int j = 0; j < igraph_matrix_ncol(K); ++j) {
      igraph_matrix_set(K, i, j, std::max<float>(Normal(mean, sd), 0));
    }
  }
}

igraph_matrix_t R_NETWORK(const char* network) {
  igraph_t net;
  FILE *ifile;
  ifile = fopen(network, "r");
  igraph_read_graph_graphml(&net, ifile, 0);
  fclose(ifile);
  igraph_matrix_t habitat;
  igraph_matrix_init(&habitat, igraph_vcount(&net), igraph_vcount(&net));
  igraph_get_adjacency(&net, &habitat, IGRAPH_GET_ADJACENCY_BOTH, 0);
  igraph_destroy(&net);  
  return habitat;
}

void SIM(igraph_matrix_t* comp_mat, igraph_matrix_t* habitat, int u, int nt, int T, float g, float r, float mr, float er, float mean, float sd, float rho, float delta, std::ostream& os) {
  igraph_matrix_t N, K;
  int N_pat = igraph_matrix_ncol(habitat);
  igraph_matrix_init(&N, N_pat, N_pat * nt * nt);
  igraph_matrix_init(&K, N_pat, N_pat * nt * nt);
  igraph_vector_t T1, T2, T3;
  igraph_vector_init(&T1, 0);
  igraph_vector_init(&T2, 0);
  igraph_vector_init(&T3, 0);
  Initializaiton(&T1, &T2, &T3, &K, &N, nt, N_pat, mean, sd);
  for(int i = 0; i < T; ++i) {
    os << u << "\t" << rho << "\t" << RICHNESS(&N, delta) << "\t" << i << "\n";
    GROWTH(comp_mat, &N, &K, &T1, &T2, r, nt);
    SPECIATION(&N, &T1, &T2, &T3, mr, nt);
    DISPERSAL(&N, rho, habitat);
    EXTINCTION(&N, er);
  }
  igraph_matrix_destroy(&N);
  igraph_matrix_destroy(&K);
  igraph_vector_destroy(&T1);
  igraph_vector_destroy(&T2);
  igraph_vector_destroy(&T3);
}

// void S_sub(int iNet, int iRep, int iRho) {
//   int delta = 10, nt = 5, Time = 15000;
//   float er = 1e-3, g = 1.5, r = 1, mr = 2.5 * 1e-5, mean = 1e3, sd = 10;
//   std::vector<float> rho = {1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
//   char FN[200] = "ahn_";
//   char buf_iNet[100];
//   sprintf(buf_iNet, "%d", iNet);
//   strcat(FN, buf_iNet);
//   char AHN[200];
//   strcpy(AHN, FN);
//   strcat(AHN, ".graphml");
//   char fn_res[200];
//   strcpy(fn_res, FN);
//   strcat(fn_res, "_rho_");
//   char buf_iRho[100];
//   sprintf(buf_iRho, "%d", iRho);
//   strcat(fn_res, buf_iRho);
//   strcat(fn_res, "_rep_");
//   char buf_iRep[100];
//   sprintf(buf_iRep, "%d", iRep);
//   strcat(fn_res, buf_iRep);
//   strcat(fn_res, "_sub.txt");
//   igraph_matrix_t habitat = R_NETWORK(AHN);
//   igraph_matrix_t sub_comp_mat = SUB_MAT(g, g, nt);
//   std::ofstream RES;
//   RES.open(fn_res);
//   RES << "rep\trho\trichness\tt\n";
//   SIM(&sub_comp_mat, &habitat, iRep, nt, Time, g, r, mr, er, mean, sd, rho[iRho], delta, RES);
//   RES.close();
//   igraph_matrix_destroy(&sub_comp_mat);
//   igraph_matrix_destroy(&habitat);
// }

void S_uns(int iNet, int iRep, int iRho) {
  int delta = 10, nt = 5, Time = 15000;
  float w_1 = .5, w_2 = .5, er = 1e-3, g = 1.5, r = 1, mr = 2.5 * 1e-5, mean = 1e3, sd = 10;
  std::vector<float> rho = {
  1e-07, 2e-07, 3e-07, 4e-07, 5e-07, 6e-07, 7e-07, 8e-07, 9e-07,
  1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06, 7e-06, 8e-06, 9e-06,
  1e-05, 2e-05, 3e-05, 4e-05, 5e-05, 6e-05, 7e-05, 8e-05, 9e-05,
  1e-04, 2e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 8e-04, 9e-04,
  1e-03
  };
  char FN[200] = "ahn_";
  char buf_iNet[100];
  sprintf(buf_iNet, "%d", iNet);
  strcat(FN, buf_iNet);
  char AHN[200];
  strcpy(AHN, FN);
  strcat(AHN, ".graphml");
  char fn_res[200];
  strcpy(fn_res, FN);
  strcat(fn_res, "_rho_");
  char buf_iRho[100];
  sprintf(buf_iRho, "%d", iRho);
  strcat(fn_res, buf_iRho);
  strcat(fn_res, "_rep_");
  char buf_iRep[100];
  sprintf(buf_iRep, "%d", iRep);
  strcat(fn_res, buf_iRep);
  strcat(fn_res, "_uns.txt");
  igraph_matrix_t habitat = R_NETWORK(AHN);
  igraph_matrix_t uns_comp_mat = UNS_MAT(g, g, nt, w_1, w_2);
  std::ofstream RES;
  RES.open(fn_res);
  RES << "rep\trho\trichness\tt\n";
  SIM(&uns_comp_mat, &habitat, iRep, nt, Time, g, r, mr, er, mean, sd, rho[iRho], delta, RES);
  RES.close();
  igraph_matrix_destroy(&uns_comp_mat);
  igraph_matrix_destroy(&habitat);
}