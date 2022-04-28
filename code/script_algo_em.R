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
GetRealParam = function(listOfSpecies){
  v_index = c()
  for (species in listOfSpecies){
    index_sp = which(nest_data$bird_names == species)
    v_index = append(v_index, index_sp)
  }
  return(nest_data[v_index,])
}

# Choix de 2 espèces European Goldfinch et Ring Ouzel afin d'avoir un mélange
# de 2 gaussiennes
df_th = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel")
                        ,proportion_alpha = c(0.3,0.7), mean_th = c(38, 298.6),
                          sd_th = c(9.1, 125.1))


# Dataframe des données d'initialisations choisies par l'utilisateur
# On a choisi les espece suivantes: "European Goldfinch"
# et "Ring Ouzel"
data_init = data.frame(alpha_init = c(0.2,0.8), mean_init = c(50, 280),
                       sd_init = c(11, 130))

simulation = function(data, n=100){
  X = NULL
  X = matrix(nrow=n,ncol=1)
  J = dim(data)[1]
  for(i in 1:n){
    Z = runif(1)
    if (Z <= data$proportion_alpha[1]){
      X[i] = rnorm(1, data$mean_th[1], data$sd_th[1])
    }else{
      k = 1
      l = 2
      Bool = FALSE
      vec_alpha = data$proportion_alpha
      cumul_alpha = cumsum(vec_alpha)
      while((Bool == FALSE) & (k < J) & (l < J+1)){
        if((cumul_alpha[k]<=Z) & (cumul_alpha[l]>=Z)){
          Bool = TRUE
          param_index = l
        }
        k = k+1
        l = l+1
      }
      X[i] = rnorm(1, data$mean_th[param_index], data$sd_th[param_index])
    }
  }
  return(X)
}

# Fonction générant le dataframe (pour l'instant vide)
# qui contiendra les données calculées à l'étape E
dataStateE = function(n, J){
  df_stateE = data.frame(matrix(NA, nrow = n, ncol = J*2+1))
  for(e in 1:dim(df_stateE)[2]){
    if(e <= J){
      names(df_stateE)[e] = paste("alpha", e, "Xi", sep="")
    }else if(e == J+1){
      names(df_stateE)[e] = "sumAlpha_j_Xi"
      index = e
    }else{
      names(df_stateE)[e] = paste("H", e-index, sep="_")
    }
  }
  return(df_stateE)
}

algo_EM = function(df, X, N=30){
  J = dim(df)[1]
  df_stateE = dataStateE(dim(X)[1], J)
  new_df = data.frame(alpha_EM = rep(NA, J), mu_EM = rep(NA, J),
                      sd_EM = rep(NA, J))
  for(n in 1:N){
    v = c()
    for(col in 1:dim(df_stateE)[2]){
      # Etape E
      colname = colnames(df_stateE)[col]
      # dans le if on calcule les alpha(j)X_i
      if(col <= J){
        df_stateE[[colname]] = df$alpha_init[col]*dnorm(X, df$mean_init[col],
                                                      df$sd_init[col])
        v = append(v, colname)
        # dans le "else if" on calcule la somme des alpha(l)X_i avec l ds [1:J]
      }else if(col == J+1){
        df_stateE[[colname]] = rowSums(df_stateE[v])
        index = col
        # dans le else on calcule les H_ij
      }else{
        ind_a = col - index
        df_stateE[[colname]] = df_stateE[[ind_a]]/(df_stateE[[J+1]])
      }
    }
    # Etape M
    for(j in 1:3){
      k = J+2
      for(i in 1:J){
        # dans le if on met à jour les alpha
        if(j == 1){
          H = df_stateE[[k]]
          new_df[i,j] = mean(H)
          k = k+1
          # dans le else if on met à jour les mu
        }else if(j == 2){
          H = df_stateE[[k]]
          new_df[i,j] = sum(H * X)/(sum(H))
          k = k+1
          # dans le else on met à jour les sigma
        }else{
          H = df_stateE[[k]]
          new_df[i,j] = sqrt( sum(H*(X-new_df$mu_EM[i])^2) / (sum(H)) )
          k = k+1
        }
      } 
    }
  }
  return(new_df)
}

# Fonction calculant la racine carré de l'erreur quadratique moyenne pour
# les  paramètre à savoir alpha, mu et sigma
# les RMSE respectives de ces 3 paramatres sont stockées dans un df
calcul_rmse = function(df_th, new_df){
  df_rmse = data.frame(rmse_alpha = NA, rmse_mu = NA, rmse_sigma = NA)
  df_th[[1]] = NULL
  for(i in 1:3){
    df_rmse[[i]] = sqrt(mean((df_th[[i]] - new_df[[i]])^2))
  }
  return(df_rmse)
}

monteCarlo = function(df_th, df_init, X, k){
  df_monteCarlo = data.frame(rmse_alpha = rep(NA, k),
                             rmse_mu = rep(NA, k),
                             rmse_sigma = rep(NA, k))
  for(i in 1:k){
    new_df = algo_EM(df_init, X, i)
    df_rmse = calcul_rmse(df_th, new_df)
    v_rmse = as.numeric(df_rmse)
    df_monteCarlo[i, ] = v_rmse
  }
  return(df_monteCarlo)
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

data_init3 = data.frame(alpha_init = c(0.4,0.5,0.1),
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

data_init4 = data.frame(alpha_init = c(0.2,0.2,0.3, 0.3),
                    mean_init = c(60, 40, 270, 313),
                    sd_init = c(7, 8.25, 42, 142))

# teste des fonctions simulation et algo_EM pour ce mélange de 4 gaussiennes
X4 = simulation(df_th4, 100)
print(df_th4)
algo_EM(data_init4, X4, 30)
