bird_names = c("European Goldfinch", "Common Linnet", "Common Chaffinch",
               "European Greenfinch", "Eurasian Bullfinch", "Hawfinch",
               "Stonechat", "European Robin", "Whinchat", "Song Thrush",
               "Common Blackbird", "Ring Ouzel", "Mistle Thrush")

mean_volume = c(38.0, 60.9, 58.3, 74.5, 45.0, 71.6, 91.0, 68.4, 51.9, 288.9,
                293.6, 298.6, 266.1)

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

# Choix de 2 espèces European Goldfinch et Ring Ouzel afin d'avoir un mélange
# de 2 gaussiennes
data_test_th = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel")
                   ,proportion_alpha = c(0.3,0.7), mean = c(38, 298.6),
                   sd = c(9.1, 125.1))


# Dataframe des données d'initialisations choisies par l'utilisateur
# On a choisi les espece suivantes: "European Goldfinch"
# et "Ring Ouzel"
data_test_init = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel"),
                       alpha_init = c(0.2,0.8), mean_init = c(50, 280),
                       sd_init = c(11, 130))

#data_th = GetRealParam
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
param_State_E = function(n, J){
  data_stateE = data.frame(matrix(NA, nrow = n, ncol = J))
  for(j in 1:J){
      names(data_stateE)[j] = paste("H(.,",j,")",sep="")
    }
  return(data_stateE)
}

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

data_test_th
X = simulation(data_test_th,100)
algo_EM(data_test_init,X,50)


# Calculs des erreurs
calcul_error = function(data_th, data_EM){
  J = dim(data_th)[1]
  df_error = data.frame(error_alpha = rep(NA, J), error_mu = rep(NA, J),
                       error_sigma = rep(NA, J))
  for(c in 2:4){
    df_error[,c-1] = abs(data_th[,c] - data_EM[,c])
  }
  return(df_error)
}


# Monte Carlo
# réalisation d'un monte carlo avec k echantillons de taille n
# avec N iterations de l'algorithme EM
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


# Choix de 2 espèces European Goldfinch et Ring Ouzel afin d'avoir un mélange
# de 2 gaussiennes
df_th2 = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel")
                    ,proportion_alpha = c(0.3,0.7), mean_th = c(38, 298.6),
                    sd_th = c(9.1, 125.1))


# Dataframe des données d'initialisations choisies par l'utilisateur
# On a choisi les espèces suivantes: "European Goldfinch"
# et "Ring Ouzel"
data_init2 = data.frame(alpha_init = c(0.2,0.8), mean_init = c(50, 280),
                        sd_init = c(11, 130))

# teste des fonctions simulation et algo_EM pour ce mélange de 2 gaussiennes
X2 = simulation(df_th2, 100)
print(df_th2)
algo_EM(data_init2, X2, 30)

# Autre exemple avec 3 espèces
df_th3 = data.frame(bird_names = c("Stonechat", "European Goldfinch",
                                   "Common Blackbird")
                    ,proportion_alpha = c(0.3,0.5,0.2),
                    mean_th = c(91.0, 38.0, 293.6),
                    sd_th = c(46.5, 9.1, 78.5))

data_init3 = data.frame(bird_names = c("Stonechat", "European Goldfinch",
                                       "Common Blackbird")
                        ,alpha_init = c(0.4,0.5,0.1),
                        mean_init = c(110, 20, 275),
                        sd_init = c(55, 21, 98.6))

# teste des fonctions simulation et algo_EM pour ce mélange de 3 gaussiennes
X3 = simulation(df_th3, 100)
print(df_th3)
algo_EM(data_init3, X3, 30)

# Autre exemple avec 4 espèces
df_th4 = data.frame(bird_names = c("Eurasian Bullfinch", "Common Chaffinch",
                                   "Mistle Thrush", "Ring Ouzel")
                    ,proportion_alpha = c(0.3,0.4,0.2, 0.1),
                    mean_th = c(45, 58.3, 266.1, 298.6),
                    sd_th = c(3.8, 15.0, 56.6, 125.1))

data_init4 = data.frame(bird_names = c("Eurasian Bullfinch", "Common Chaffinch",
                                       "Mistle Thrush", "Ring Ouzel"),
                        alpha_init = c(0.2,0.2,0.3, 0.3),
                        mean_init = c(60, 40, 270, 313),
                        sd_init = c(7, 8.25, 42, 142))

# teste des fonctions simulation et algo_EM pour ce mélange de 4 gaussiennes
X4 = simulation(df_th4, 100)
print(df_th4)
algo_EM(data_init4, X4, 30)




# Boxplots

vect1 <- seq(1, 298, by = 3) # vecteur 1
vect2 <- seq(2, 299, by = 3) # vecteur 2
vect3 <- seq(3, 300, by = 3) # vecteur 3



# Pour les \alpha
df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
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
df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
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
df_monteCarlo_50 = Monte_Carlo(df_th3, data_init3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(df_th3, data_init3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(df_th3, data_init3, 100, 500, 10)
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



# Pour les \v_j
df_monteCarlo_50 = Monte_Carlo(data_test_th, data_test_init, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(data_test_th, data_test_init, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(data_test_th, data_test_init, 100, 500, 10)
par(mfrow = c(2,3))
boxplot(df_monteCarlo_50[vectp, 4]) # \v_1
boxplot(df_monteCarlo_100[vectp, 4]) # \v_1
boxplot(df_monteCarlo_500[vectp, 4]) # \v_1
boxplot(df_monteCarlo_50[vecti, 4]) # \v_2
boxplot(df_monteCarlo_100[vecti, 4]) # \v_2
boxplot(df_monteCarlo_500[vecti, 4]) # \v_2


