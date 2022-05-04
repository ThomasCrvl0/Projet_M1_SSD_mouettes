# Vecteur dans lequel on a stocké les noms des 13 espèces d'oiseaux présentent dans notre document
bird_names = c("European Goldfinch", "Common Linnet", "Common Chaffinch",
               "European Greenfinch", "Eurasian Bullfinch", "Hawfinch",
               "Stonechat", "European Robin", "Whinchat", "Song Thrush",
               "Common Blackbird", "Ring Ouzel", "Mistle Thrush")

# Vecteur dans lequel on a stocké les moyennes respectives des volumes des nids des 13 espèces d'oiseaux
mean_volume = c(38.0, 60.9, 58.3, 74.5, 45.0, 71.6, 91.0, 68.4, 51.9, 288.9,
                293.6, 298.6, 266.1)
# Vecteur dans lequel on a stocké les écarts-types respectifs des volumes des nids des 13 espèces d'oiseaux
sd_volume = c(9.1,  20.8, 15.0,  12.2,  3.8,  12.9,  46.5, 29.8, 27.4, 55.9,
              78.5,  125.1,  56.6)

# Création du dataframe contenant les données sur les 13 nids d'oiseaux
# Chaque ligne du df correspond à une espèce d'oiseau
# Chaque espèce est associée à 2 colonnes qui nous renseignent sur
# le volume moyen des nids associée à cette espèce respectivement
# l'ecart type des nids associées à cette espèce
nest_data = data.frame(bird_names, mean_volume, sd_volume)

# On crée le dataframe regroupant les moyennes et les écarts-types théoriques
# associés aux espèces d'oiseaux choisies par l'utilisateur
# Les espèces choisies sont stockées dans un vecteur pris en paramètre
GetRealParam = function(list_Species){
  v_index = c()
  for (species in list_Species){
    index_sp = which(nest_data$bird_names == species)
    v_index = append(v_index, index_sp)
  }
  return(nest_data[v_index,])
}


# Cette fonction génère aléatoirement un échantillon de taille n issue d'un mélange gaussien
# Les paramètres (alpha, mu et sigma) des différents mélanges gaussiens sont contenus dans data_th
# La fonction prends en argument data_th (le dataframe contenant les paramètres alpha, mu et sigma)
# et n qui indique le nombre de valeurs que l'on souhaite générer aléatoirement
# La fonction retourne un vecteur contenant les n valeurs du mélange gaussien
# qui ont été generées aléatoirement
simulation = function(data_th, n=100){
  X = rep(NA,n) #echantillon
  J = dim(data_th)[1] #nb de mélange
  vect_alpha = data_th[,2]
  vect_mean = data_th[,3]
  vect_sd = data_th[,4]
  for(i in 1:n){
    Z = runif(1)
    if (Z <= vect_alpha[1]){
      X[i] = rnorm(1, vect_mean[1], vect_sd[1])
    }else{
      k = 1
      l = 2
      Bool = FALSE
      cumul_alpha = cumsum(vect_alpha)
      while(Bool == FALSE){
        if((cumul_alpha[k]<=Z) & (cumul_alpha[l]>=Z)){
          Bool = TRUE
          param_index = l
        }
        k = k+1
        l = l+1
      }
      X[i] = rnorm(1, vect_mean[param_index], vect_sd[param_index])
    }
  }
  return(X)
}

# Fonction générant dataframe (vide)
# contenant les valeurs de P_thetat(j|X = X_i) (etape E de l'algo)
# Cette fonction prend en argument n (la taille de l'échantillon)
# Et J le nombre de mélange
param_State_E = function(n, J){
  data_stateE = data.frame(matrix(NA, nrow = n, ncol = J))
  for(j in 1:J){
    names(data_stateE)[j] = paste("H(.,",j,")",sep="")
  }
  return(data_stateE)
}

# Cette fonction est une implémentation de l'algorithme EM
# Elle prend en argument data_init
# (data_init est un dataframe contenant les paramètres (alpha, mu, sigma) initiaux choisis)
# un vecteur X qui est l'échantillon de taille n
# K le nombre d'itérations de l'algorithme EM
# Cette fonction va retourner un dataframe contenant les valeurs des paramètres alpha, mu et sigma
# qui ont été mises à jour après K itérations de l'algorithme EM

algo_EM = function(data_init, X, K){
  J = dim(data_init)[1]
  n = length(X)
  data_stateE = param_State_E(n, J)
  
  for(k in 1:K){
    vect_alpha = data_init[,2] #de longueur J
    vect_mean = data_init[,3]
    vect_sd = data_init[,4]
    v = rep(0,n) #vecteur contenant la somme des des numérateurs de P_thetat(j|X = X_i) 
    #pour chaque valeur de l'echantillon 
    
    # Etape E
    for(j in 1:J){
      # On remplie le taleau param_State_E contenant P_thetat(j|X=X_i)
      data_stateE[,j] = vect_alpha[j]*dnorm(X,vect_mean[j],vect_sd[j])
      v = v+data_stateE[,j]
    }
    for(j in 1:J){
      data_stateE[,j] = data_stateE[,j]/v
    }
    
    # Etape M
    H = data_stateE
    for(col in 2:4){ #on met a jour le data_init
      for(ind in 1:J){
        
        # on met à jour les alpha
        if(col == 2){
          data_init[,col][ind] = mean(H[,ind])
        }
        # on met à jour les mu
        if(col == 3){
          data_init[,col][ind] = (sum(X*H[,ind]))/(sum(H[,ind]))
        } 
        # on met à jour les sigma
        if(col == 4){
          data_init[,col][ind] = sqrt((sum( (X-rep(data_init[,col-1][ind],n))^2
                                            *H[,ind] ))/sum(H[,ind]))
        }
      }
    }
  }
  new_df = data_init
  return(new_df)
}


# Cette fonction calcule l'erreur absolue de tous les paramètres alpha, mu et sigma
# elle calcule la différence en valeur absolue entre la valeur du paramètre estimé
# par l'algorithme EM et la valeur du paramètre théorique
# Elle prend en paramètre data_th le dataframe contenant les paramètres alpha, mu et sigma théoriques
# et Data_EM dataframe contenant les paramètres alpha, mu et sigma estimés par l'algorithme EM
# Cette fonction retourne le dataframe contenant les erreurs absolues de chaque paramètre

calcul_error = function(data_th, data_EM){
  J = dim(data_th)[1]
  df_error = data.frame(error_alpha = rep(NA, J), error_mu = rep(NA, J),
                        error_sigma = rep(NA, J))
  for(c in 2:4){
    df_error[,c-1] = abs(data_th[,c] - data_EM[,c])
  }
  return(df_error)
}


# Cette fonction calcule l'erreur absolue entre les paramètres estimés par l'algo EM
# et les paramètres théoriques pour k echantillons générés aléatoirement
# Elle prend en argument data_th le dataframe contenant les paramètres alpha, mu et sigma théoriques
# data_init le dataframe contenant les paramètres (alpha, mu, sigma) initiaux choisis
# k le nombre d'échantillons que l'on souhaite générer aléatoirement
# n la taille des échantillons (les k échantillons seront de tailles n)
# N le nombre d'iterations de l'algorithme EM
# Elle retourne un dataframe contenant les erreurs absolues de chaque paramètre
# des k échantillons de Monte Carlo
Monte_Carlo = function(data_th, data_init, k, n, N){
  J = dim(data_th)[1]
  df_MonteCarlo = data.frame()
  iteration = c()
  for(i in 1:k){
    X = simulation(data_th, n)
    data_EM = algo_EM(data_init, X, N)
    df_error = calcul_error(data_th, data_EM)
    df_MonteCarlo = rbind((df_MonteCarlo), df_error)
    v_iter = rep(paste("itération",i,sep="_"), J)
    iteration = append(iteration, v_iter)
  }
  df_MonteCarlo = cbind(iteration, df_MonteCarlo)
  return(df_MonteCarlo)
}

#############################################################################
# partie tests
#############################################################################
# On teste l'algorithme EM sur un mélange de 2 gaussiennes

# Choix de 2 espèces European Goldfinch et Ring Ouzel afin d'avoir un mélange
# de 2 gaussiennes
data_test_th = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel")
                          ,proportion_alpha = c(0.3,0.7), mean = c(38, 298.6),
                          sd = c(9.1, 125.1))


# Dataframe des données d'initialisations choisies par l'utilisateur
# On a choisi les espèces suivantes: "European Goldfinch"
# et "Ring Ouzel"
data_test_init = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel"),
                            alpha_init = c(0.2,0.8), mean_init = c(50, 280),
                            sd_init = c(11, 130))

data_test_th
X = simulation(data_test_th,100)
algo_EM(data_test_init,X,50)


# Exemple avec 3 espèces ayant des variables bien séparées
good_dfTestTh3 = data.frame(bird_names = c("European Goldfinch", "Stonechat", "Ring Ouzel")
                    ,proportion_alpha = c(0.3,0.5,0.2),
                    mean_th = c(38.0, 91.0, 298.6),
                    sd_th = c(9.1, 46.5, 125.1))

good_dfTestInit3 = data.frame(bird_names = c("European Goldfinch", "Stonechat", "Ring Ouzel")
                              ,proportion_alpha = c(0.35,0.57,0.13),
                              mean_th = c(44.4, 100, 278.3),
                              sd_th = c(10, 57, 130))


# teste des fonctions simulation et algo_EM pour ce mélange de 3 gaussiennes avec forte séparation
X3 = simulation(good_dfTestTh3, 100)
print(good_dfTestTh3)
algo_EM(good_dfTestInit3, X3, 30)

# on affiche les boxplots pour ce melange de " gausiennes avec forte séparation
# Boxplots

vect1 <- seq(1, 298, by = 3) # vecteur 1
vect2 <- seq(2, 299, by = 3) # vecteur 2
vect3 <- seq(3, 300, by = 3) # vecteur 3



# Pour les \alpha
df_monteCarlo_50 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_100[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_500[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_50[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_100[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_500[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_50[vect3, 2]) # \alpha_3
boxplot(df_monteCarlo_100[vect3, 2]) # \alpha_3
boxplot(df_monteCarlo_500[vect3, 2]) # \alpha_3



# Pour les \mu
df_monteCarlo_50 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_100[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_500[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_50[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_100[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_500[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_50[vect3, 3]) # \mu_3
boxplot(df_monteCarlo_100[vect3, 3]) # \mu_3
boxplot(df_monteCarlo_500[vect3, 3]) # \mu_3

# Pour les \sigma
df_monteCarlo_50 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_100[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_500[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_50[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_100[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_500[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_50[vect3, 4]) # \sigma_3
boxplot(df_monteCarlo_100[vect3, 4]) # \sigma_3
boxplot(df_monteCarlo_500[vect3, 4]) # \sigma_3

# Exemple avec 3 espèces ayant des variables mal séparées
bad_dfTestTh3 = data.frame(bird_names = c("European Greenfinch", "Hawfinch", "Common Chaffinch")
                            ,proportion_alpha = c(0.28,0.12,0.6),
                            mean_th = c(74.5, 71.6, 58.3),
                            sd_th = c(12.2, 12.9, 15.0))

bad_dfTestInit3 = data.frame(bird_names = c("European Greenfinch", "Hawfinch", "Common Chaffinch")
                              ,proportion_alpha = c(0.3,0.2,0.5),
                              mean_th = c(78.3, 69.7, 60),
                              sd_th = c(11.8, 13.3, 17))


# teste des fonctions simulation et algo_EM pour ce mélange de 3 gaussiennes avec faible séparation
X3_bad = simulation(bad_dfTestTh3, 100)
print(bad_dfTestTh3)
algo_EM(bad_dfTestInit3, X3_bad, 30)

# on affiche les boxplots pour ce melange de gausiennes avec faible séparation
# Boxplots

vect1 <- seq(1, 298, by = 3) # vecteur 1
vect2 <- seq(2, 299, by = 3) # vecteur 2
vect3 <- seq(3, 300, by = 3) # vecteur 3



# Pour les \alpha
df_monteCarlo_50 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_100[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_500[vect1, 2]) # \alpha_1
boxplot(df_monteCarlo_50[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_100[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_500[vect2, 2]) # \alpha_2
boxplot(df_monteCarlo_50[vect3, 2]) # \alpha_3
boxplot(df_monteCarlo_100[vect3, 2]) # \alpha_3
boxplot(df_monteCarlo_500[vect3, 2]) # \alpha_3



# Pour les \mu
df_monteCarlo_50 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_100[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_500[vect1, 3]) # \mu_1
boxplot(df_monteCarlo_50[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_100[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_500[vect2, 3]) # \mu_2
boxplot(df_monteCarlo_50[vect3, 3]) # \mu_3
boxplot(df_monteCarlo_100[vect3, 3]) # \mu_3
boxplot(df_monteCarlo_500[vect3, 3]) # \mu_3

# Pour les \sigma
df_monteCarlo_50 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 500, 10)
par(mfrow = c(3,3))
boxplot(df_monteCarlo_50[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_100[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_500[vect1, 4]) # \sigma_1
boxplot(df_monteCarlo_50[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_100[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_500[vect2, 4]) # \sigma_2
boxplot(df_monteCarlo_50[vect3, 4]) # \sigma_3
boxplot(df_monteCarlo_100[vect3, 4]) # \sigma_3
boxplot(df_monteCarlo_500[vect3, 4]) # \sigma_3

# teste des fonctions simulation et algo_EM pour ce mélange de 4 gaussiennes
X4 = simulation(df_th4, 100)
print(df_th4)
algo_EM(data_init4, X4, 30)




# Boxplots

# vect1 <- seq(1, 298, by = 3) # vecteur 1
# vect2 <- seq(2, 299, by = 3) # vecteur 2
# vect3 <- seq(3, 300, by = 3) # vecteur 3



# Pour les \alpha
# df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
# df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
# df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
# par(mfrow = c(3,3))
# boxplot(df_monteCarlo_50[vect1, 2]) # \alpha_1
# boxplot(df_monteCarlo_100[vect1, 2]) # \alpha_1
# boxplot(df_monteCarlo_500[vect1, 2]) # \alpha_1
# boxplot(df_monteCarlo_50[vect2, 2]) # \alpha_2
# boxplot(df_monteCarlo_100[vect2, 2]) # \alpha_2
# boxplot(df_monteCarlo_500[vect2, 2]) # \alpha_2
# boxplot(df_monteCarlo_50[vect3, 2]) # \alpha_3
# boxplot(df_monteCarlo_100[vect3, 2]) # \alpha_3
# boxplot(df_monteCarlo_500[vect3, 2]) # \alpha_3



# Pour les \mu
# df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
# df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
# df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
# par(mfrow = c(3,3))
# boxplot(df_monteCarlo_50[vect1, 3]) # \mu_1
# boxplot(df_monteCarlo_100[vect1, 3]) # \mu_1
# boxplot(df_monteCarlo_500[vect1, 3]) # \mu_1
# boxplot(df_monteCarlo_50[vect2, 3]) # \mu_2
# boxplot(df_monteCarlo_100[vect2, 3]) # \mu_2
# boxplot(df_monteCarlo_500[vect2, 3]) # \mu_2
# boxplot(df_monteCarlo_50[vect3, 3]) # \mu_3
# boxplot(df_monteCarlo_100[vect3, 3]) # \mu_3
# boxplot(df_monteCarlo_500[vect3, 3]) # \mu_3

# Pour les \sigma
# df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
# df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
# df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
# par(mfrow = c(3,3))
# boxplot(df_monteCarlo_50[vect1, 4]) # \sigma_1
# boxplot(df_monteCarlo_100[vect1, 4]) # \sigma_1
# boxplot(df_monteCarlo_500[vect1, 4]) # \sigma_1
# boxplot(df_monteCarlo_50[vect2, 4]) # \sigma_2
# boxplot(df_monteCarlo_100[vect2, 4]) # \sigma_2
# boxplot(df_monteCarlo_500[vect2, 4]) # \sigma_2
# boxplot(df_monteCarlo_50[vect3, 4]) # \sigma_3
# boxplot(df_monteCarlo_100[vect3, 4]) # \sigma_3
# boxplot(df_monteCarlo_500[vect3, 4]) # \sigma_3



# Pour les \v_j
# df_monteCarlo_50 = Monte_Carlo(data_test_th, data_test_init, 100, 50, 10)
# df_monteCarlo_100 = Monte_Carlo(data_test_th, data_test_init, 100, 100, 10)
# df_monteCarlo_500 = Monte_Carlo(data_test_th, data_test_init, 100, 500, 10)
# par(mfrow = c(2,3))
# boxplot(df_monteCarlo_50[vectp, 4]) # \v_1
# boxplot(df_monteCarlo_100[vectp, 4]) # \v_1
# boxplot(df_monteCarlo_500[vectp, 4]) # \v_1
# boxplot(df_monteCarlo_50[vecti, 4]) # \v_2
# boxplot(df_monteCarlo_100[vecti, 4]) # \v_2
# boxplot(df_monteCarlo_500[vecti, 4]) # \v_2
