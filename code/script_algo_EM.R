library(ggplot2)
library(ggpubr)

# Vecteur dans lequel on a stocké les noms des 13 espèces d'oiseaux présentent 
# dans notre document
bird_names = c("European Goldfinch", "Common Linnet", "Common Chaffinch",
               "European Greenfinch", "Eurasian Bullfinch", "Hawfinch",
               "Stonechat", "European Robin", "Whinchat", "Song Thrush",
               "Common Blackbird", "Ring Ouzel", "Mistle Thrush")

# Vecteur dans lequel on a stocké les moyennes respectives des volumes des nids
# des 13 espèces d'oiseaux
mean_volume = c(38.0, 60.9, 58.3, 74.5, 45.0, 71.6, 91.0, 68.4, 51.9, 288.9,
                293.6, 298.6, 266.1)
# Vecteur dans lequel on a stocké les écarts-types respectifs des volumes des
# nids des 13 espèces d'oiseaux
sd_volume = c(9.1,  20.8, 15.0,  12.2,  3.8,  12.9,  46.5, 29.8, 27.4, 55.9,
              78.5,  125.1,  56.6)

# Création du dataframe contenant les données sur les 13 nids d'oiseaux
# Chaque ligne du df correspond à une espèce d'oiseau
# Chaque espèce est associée à 2 colonnes qui nous renseignent sur
# le volume moyen des nids associée à cette espèce respectivement
# l'ecart type des nids associées à cette espèce
nest_data = data.frame(bird_names, mean_volume, sd_volume)

# Cette fonction prend en argument df le dataframe regroupant toutes les données
# c'est à dire les noms, les moyennes et les écart-types des nids des 13
# espèces d'oiseaux
# n le nombre d'espèces à tirer aléatoirement
# cette fonction retourne un dataframe stockant les noms,
# moyennes et écart-types des nids des n espèces d'oiseaux tirées aléatoirement
random_species = function(df, n){
  J = dim(df)[1]
  index_species = sample(1:J, n, replace = FALSE)
  random_prop = runif(n, 0, 1)
  random_prop = random_prop/sum(random_prop)
  random_df = df[index_species,]
  random_df = cbind(proportion_alpha = random_prop, random_df)
  random_df = random_df[,c(2,1,3,4)]
  return(random_df)
}

# Cette fonction affiche graphiquement la distribution d'un mélange gaussien
# Elle prend en argument df le dataframe contenant les données 
# (noms, moyennes, ecarts-types) des nids des espèces d'oiseaux 
# et X l’échantillon du mélange gaussien
plot_distrib = function(df, X){
  data_distrib = data.frame(gaussian_mixture = X)
  ggplot(data_distrib, aes(x=data_distrib[,'gaussian_mixture'])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 100) +
  geom_density(alpha=.3, color = "aquamarine2", fill="red", size=1.5) + 
  ggtitle("Distribution d'un mélange gaussien")
}

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


# Cette fonction génère aléatoirement un échantillon de taille n issue 
# d'un mélange gaussien
# Les paramètres (alpha, mu et sigma) des différents mélanges gaussiens 
# sont contenus dans data_th
# La fonction prends en argument data_th 
# (le dataframe contenant les paramètres alpha, mu et sigma)
# et n qui indique le nombre de valeurs que l'on souhaite générer aléatoirement
# La fonction retourne un vecteur contenant les n valeurs du mélange gaussien
# qui ont été generées aléatoirement
simulation = function(data_th, n=100){
  X = rep(NA,n) #echantillon
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
# (data_init est un dataframe contenant les paramètres (alpha, mu, sigma)
# initiaux choisis)
# un vecteur X qui est l'échantillon de taille n
# K le nombre d'itérations de l'algorithme EM
# Cette fonction va retourner un dataframe contenant les valeurs des
# paramètres alpha, mu et sigma
# qui ont été mises à jour après K itérations de l'algorithme EM

algo_EM = function(data_init, X, K){
  J = dim(data_init)[1]
  n = length(X)
  data_stateE = param_State_E(n, J)
  
  for(k in 1:K){
    vect_alpha = data_init[,2] #de longueur J
    vect_mean = data_init[,3]
    vect_sd = data_init[,4]
    # vecteur contenant la somme des des numérateurs de P_thetat(j|X = X_i)
    # pour chaque valeur de l’échantillon 
    v = rep(0,n) 
    
    # Etape E
    for(j in 1:J){
      # On remplie le tableau param_State_E contenant P_thetat(j|X=X_i)
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
  colnames(new_df) = c('bird_names', 'alpha', 'mu', 'sigma')
  return(new_df)
}


# Cette fonction calcule l'erreur absolue de tous les paramètres alpha, mu et sigma
# elle calcule la différence en valeur absolue entre la valeur du paramètre estimé
# par l'algorithme EM et la valeur du paramètre théorique
# Elle prend en paramètre data_th le dataframe contenant les 
# paramètres alpha, mu et sigma théoriques
# et Data_EM dataframe contenant les paramètres alpha, mu et sigma estimés 
# par l'algorithme EM
# Cette fonction retourne le dataframe contenant les erreurs absolues de chaque
# paramètre

calcul_error = function(data_th, data_EM){
  J = dim(data_th)[1]
  df_error = data.frame(error_alpha = rep(NA, J), error_mu = rep(NA, J),
                        error_sigma = rep(NA, J))
  for(c in 2:4){
    df_error[,c-1] = abs(data_th[,c] - data_EM[,c])
  }
  return(df_error)
}


# Cette fonction calcule l'erreur absolue entre les paramètres estimés par
# l'algo EM
# et les paramètres théoriques pour k echantillons générés aléatoirement
# Elle prend en argument data_th le dataframe contenant les paramètres
# alpha, mu et sigma théoriques
# data_init le dataframe contenant les paramètres (alpha, mu, sigma)
# initiaux choisis
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


#Paramètres initiaux 

param_quantile1 <- function(X,J){
  N <- length(X)
  Alpha = rep(1/J,J)
  Mean <- rep(NA,J)
  
  for(j in 1:J){
    Mean[j] <- floor(j*N/J+1)
  }
  
  X <- sort(X)
  Sd <- sqrt(var(X[1:Mean[1]]))
  Sd <- rep(Sd,J)
  
  data <- data.frame(init_alpha = Alpha, init_mu = Mean, init_sd = Sd)
  return(data)
}

param_quantile2 <- function(X,J){
  N <- length(X)
  Alpha = rep(1/J,J)
  Mean <- rep(NA,J)
  X <- sort(X)
  
  index1_Q <- 1
  index2_Q <- 0
  for(j in 1:J-1){
    index2_Q <- floor(j*N/J)+1
    Mean[j] <- mean(X[index1_Q:index2_Q])
    index1_Q <- index2_Q+1
  }
  Mean[J] <- mean(X[index1_Q:N])
  
  Sd <- sqrt(var(X[1:floor(N/J)+1]))
  Sd <- rep(Sd,J)
  
  data <- data.frame(init_alpha = Alpha, init_mu = Mean, init_sd = Sd)
  return(data)
}

param_kmeans <- function(X,J){
  density <- density(X)
  abscisse <- density$x
  ordonné <- density$y
  
  m_X <- abscisse[ggpmisc:::find_peaks(ordonné)]#positions des max
  
  centers <- data.frame(moyennes = m_X)
  size_m_X <- length(m_X)
  alpha <- rep(NA,J)
  var <- rep(NA,J)
  
  if (size_m_X > J){
    min <- m_X[m_X < mean(X)] 
    max <- m_X[m_X > mean(X)]
    min_max <- sort(c(abs(min-mean(X)),abs(max-mean(X))))
    i <- 1 
    while (size_m_X > J){
      if (!is.null(min_max)){
        m_X <- m_X[m_X != min_max[i]] #On enlève les max les plus éloignés
        size_m_X <- length(m_X)
      }
      else{
        m_X2 <- sort(m_X)
        m_X <- m_X[m_X != m_X2[i]] #On enlève les max les plus faibles
      }
      i <- i+1
    }
  }
  
  if (size_m_X==J){
    k <- kmeans(X,centers)
    m_X <- k$centers #permet d'éviter d'avoir de faux max
    indiv <- k$cluster
    for (j in 1:J){
      alpha[j] <- (1/100)*k$size[j]
      var[j] <- sqrt(var(X[indiv==j]))
    }
    param_init <- data.frame(espèces = 1:J, alpha = alpha, moyenne = m_X, variance = sqrt(var))
    return(param_init)
  }
  
  
  while (size_m_X<J){
    k <- kmeans(X,centers)
    indice_kmeans <- which.max(k$size)
    nv_kmeans <- kmeans(X[k$cluster==indice_kmeans],2)
    moyennes <- k$centers[-(indice_kmeans:size_m_X)]
    moyennes <- c(moyennes,nv_kmeans$centers)
    moyennes <- c(moyennes,k$centers[-(1:indice_kmeans)])
    m_X <- moyennes
    centers <- data.frame(moyennes = m_X)
    size_m_X <- size_m_X+1
  }
  
  k <- kmeans(X,centers)
  for (j in 1:J){
    indiv <- k$cluster
    alpha[j] <- (1/100)*k$size[j]
    var[j] <- sqrt(var(X[indiv==j]))
  }
  m_X <- k$centers
  param_init <- data.frame(espèces = 1:J, alpha = alpha, moyenne = m_X, 
                           variance = sqrt(var))
  return(param_init)
}


log_Vrais_X <- function(data_param,X){
  J <- length(data_param[,1])
  logVrai_X <- 0
  for(i in 1:length(X)){
    sum_j <- 0
    for(j in 1:J){
      sum_j <- sum_j +  data_param[j,2]*dnorm(X[i],mean = data_param[j,3], 
                                              sd = data_param[j,4])
    }
    logVrai_X <- logVrai_X + log(sum_j)
  }
  return(logVrai_X)
}


param_init <- function(data,X){
  J <- length(data[,1])
  col <- length(data[1,])-1
  res <- -Inf
  data_res <- data.frame()
  for(i in seq(1,col,by = 3)){
    data_i <- data.frame(espèces = data[,1], alpha = data[,i+1],
                         moyennes = data[,i+2], sd = data[,i+3])
    log_Vrai <- log_Vrais_X(data_i,X)
    if (res <= log_Vrai){
      res <- log_Vrai
      data_res <- data_i
    }
  }
  return(data_res)
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
good_dfTestTh3 = data.frame(bird_names = c("European Goldfinch", "Stonechat",
                                           "Ring Ouzel")
                            ,proportion_alpha = c(0.3,0.5,0.2),
                            mean_th = c(38.0, 91.0, 135),
                            sd_th = c(9.1, 46.5, 125.1))

good_dfTestInit3 = data.frame(bird_names = c("European Goldfinch", "Stonechat",
                                             "Ring Ouzel")
                              ,proportion_alpha = c(0.35,0.57,0.13),
                              mean_th = c(44.4, 100, 278.3),
                              sd_th = c(10, 57, 130))


# teste des fonctions simulation et algo_EM pour ce mélange de 3 gaussiennes 
# avec forte séparation
X3 = simulation(good_dfTestTh3, 500)
print(good_dfTestTh3)
algo_EM(good_dfTestInit3, X3, 30)

# On affiche l'histogramme et la distribution pour le mélange avec
# forte séparation
dataGoodSample = data.frame(gaussian_mixture = X3)
ggplot(dataGoodSample, aes(x=dataGoodSample[,'gaussian_mixture'])) + 
geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 100) +
geom_density(alpha=.3, color = "aquamarine2", fill="red", size=1.5) + 
ggtitle("Distribution d'un mélange gaussien à forte séparation")

# on affiche les boxplots pour ce melange de gausiennes avec forte séparation

vect1 <- seq(1, 298, by = 3) # vecteur 1
vect2 <- seq(2, 299, by = 3) # vecteur 2
vect3 <- seq(3, 300, by = 3) # vecteur 3

df_monteCarlo_50 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 50, 10)
df_monteCarlo_100 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 100, 10)
df_monteCarlo_500 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 500, 10)

# Pour les \alpha

a1 = ggplot(df_monteCarlo_50[vect1,], aes(y = df_monteCarlo_50[vect1, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="brown2")
a2 = ggplot(df_monteCarlo_100[vect1,], aes(y = df_monteCarlo_100[vect1, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="orange")
a3 = ggplot(df_monteCarlo_500[vect1,], aes(y = df_monteCarlo_500[vect1, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="aquamarine")
a4 = ggplot(df_monteCarlo_50[vect2,], aes(y = df_monteCarlo_50[vect2, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="brown2")
a5 = ggplot(df_monteCarlo_100[vect2,], aes(y = df_monteCarlo_100[vect2, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="orange")
a6 = ggplot(df_monteCarlo_500[vect2,], aes(y = df_monteCarlo_500[vect2, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="aquamarine")
a7 = ggplot(df_monteCarlo_50[vect3,], aes(y = df_monteCarlo_50[vect3, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="brown2")
a8 = ggplot(df_monteCarlo_100[vect3,], aes(y = df_monteCarlo_100[vect3, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="orange")
a9 = ggplot(df_monteCarlo_500[vect3,], aes(y = df_monteCarlo_500[vect3, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="aquamarine")
ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, ncol = 3, nrow = 3)
# save with width = 800 and height = 1200


# Pour les \mu

m1 = ggplot(df_monteCarlo_50[vect1,], aes(y = df_monteCarlo_50[vect1, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="brown2")
m2 = ggplot(df_monteCarlo_100[vect1,], aes(y = df_monteCarlo_100[vect1, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="orange")
m3 = ggplot(df_monteCarlo_500[vect1,], aes(y = df_monteCarlo_500[vect1, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="aquamarine")
m4 = ggplot(df_monteCarlo_50[vect2,], aes(y = df_monteCarlo_50[vect2, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="brown2")
m5 = ggplot(df_monteCarlo_100[vect2,], aes(y = df_monteCarlo_100[vect2, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="orange")
m6 = ggplot(df_monteCarlo_500[vect2,], aes(y = df_monteCarlo_500[vect2, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="aquamarine")
m7 = ggplot(df_monteCarlo_50[vect3,], aes(y = df_monteCarlo_50[vect3, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="brown2")
m8 = ggplot(df_monteCarlo_100[vect3,], aes(y = df_monteCarlo_100[vect3, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="orange")
m9 = ggplot(df_monteCarlo_500[vect3,], aes(y = df_monteCarlo_500[vect3, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="aquamarine")
ggarrange(m1, m2, m3, m4, m5, m6, m7, m8, m9, ncol = 3, nrow = 3)
# save with width = 800 and height = 1200

# Pour les \sigma

s1 = ggplot(df_monteCarlo_50[vect1,], aes(y = df_monteCarlo_50[vect1, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="brown2")
s2 = ggplot(df_monteCarlo_100[vect1,], aes(y = df_monteCarlo_100[vect1, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="orange")
s3 = ggplot(df_monteCarlo_500[vect1,], aes(y = df_monteCarlo_500[vect1, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="aquamarine")
s4 = ggplot(df_monteCarlo_50[vect2,], aes(y = df_monteCarlo_50[vect2, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="brown2")
s5 = ggplot(df_monteCarlo_100[vect2,], aes(y = df_monteCarlo_100[vect2, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="orange")
s6 = ggplot(df_monteCarlo_500[vect2,], aes(y = df_monteCarlo_500[vect2, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="aquamarine")
s7 = ggplot(df_monteCarlo_50[vect3,], aes(y = df_monteCarlo_50[vect3, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="brown2")
s8 = ggplot(df_monteCarlo_100[vect3,], aes(y = df_monteCarlo_100[vect3, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="orange")
s9 = ggplot(df_monteCarlo_500[vect3,], aes(y = df_monteCarlo_500[vect3, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="aquamarine")
ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, s9, ncol = 3, nrow = 3)


# Exemple avec 3 espèces ayant des variables mal séparées
bad_dfTestTh3 = data.frame(bird_names = c("European Greenfinch", "Hawfinch",
                                          "Common Chaffinch")
                           ,proportion_alpha = c(0.28,0.12,0.6),
                           mean_th = c(74.5, 71.6, 58.3),
                           sd_th = c(12.2, 12.9, 15.0))

bad_dfTestInit3 = data.frame(bird_names = c("European Greenfinch", "Hawfinch",
                                            "Common Chaffinch")
                             ,proportion_alpha = c(0.3,0.2,0.5),
                             mean_th = c(78.3, 69.7, 60),
                             sd_th = c(11.8, 13.3, 17))


# teste des fonctions simulation et algo_EM pour ce mélange de 3 gaussiennes avec faible séparation
X3_bad = simulation(bad_dfTestTh3, 100)
print(bad_dfTestTh3)
algo_EM(bad_dfTestInit3, X3_bad, 30)

# On affiche l'histogramme et la distribution pour le mélange avec
# faible séparation
dataBadSample = data.frame(gaussian_mixture = X3_bad)
ggplot(dataBadSample, aes(x=dataBadSample[,'gaussian_mixture'])) + 
geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
geom_density(alpha=.3, color = "purple", fill="deepskyblue", size=1.5) +
ggtitle("Distribution d'un mélange gaussien à faible séparation")

# on affiche les boxplots pour ce melange de gausiennes avec faible séparation

vect1 <- seq(1, 298, by = 3) # vecteur 1
vect2 <- seq(2, 299, by = 3) # vecteur 2
vect3 <- seq(3, 300, by = 3) # vecteur 3

df_monteCarloBad_50 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 50, 10)
df_monteCarloBad_100 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 100, 10)
df_monteCarloBad_500 = Monte_Carlo(bad_dfTestTh3, bad_dfTestInit3, 100, 500, 10)

# Pour les \alpha

ab1 = ggplot(df_monteCarloBad_50[vect1,], aes(y = df_monteCarloBad_50[vect1, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="brown2")
ab2 = ggplot(df_monteCarloBad_100[vect1,], aes(y = df_monteCarloBad_100[vect1, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="orange")
ab3 = ggplot(df_monteCarloBad_500[vect1,], aes(y = df_monteCarloBad_500[vect1, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha1") + geom_boxplot(fill="aquamarine")
ab4 = ggplot(df_monteCarloBad_50[vect2,], aes(y = df_monteCarloBad_50[vect2, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="brown2")
ab5 = ggplot(df_monteCarloBad_100[vect2,], aes(y = df_monteCarloBad_100[vect2, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="orange")
ab6 = ggplot(df_monteCarloBad_500[vect2,], aes(y = df_monteCarloBad_500[vect2, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha2") + geom_boxplot(fill="aquamarine")
ab7 = ggplot(df_monteCarloBad_50[vect3,], aes(y = df_monteCarloBad_50[vect3, 2])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="brown2")
ab8 = ggplot(df_monteCarloBad_100[vect3,], aes(y = df_monteCarloBad_100[vect3, 2])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="orange")
ab9 = ggplot(df_monteCarloBad_500[vect3,], aes(y = df_monteCarloBad_500[vect3, 2])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour alpha3") + geom_boxplot(fill="aquamarine")
ggarrange(ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8, ab9, ncol = 3, nrow = 3)


# Pour les \mu

mb1 = ggplot(df_monteCarloBad_50[vect1,], aes(y = df_monteCarloBad_50[vect1, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="brown2")
mb2 = ggplot(df_monteCarloBad_100[vect1,], aes(y = df_monteCarloBad_100[vect1, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="orange")
mb3 = ggplot(df_monteCarloBad_500[vect1,], aes(y = df_monteCarloBad_500[vect1, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu1") + geom_boxplot(fill="aquamarine")
mb4 = ggplot(df_monteCarloBad_50[vect2,], aes(y = df_monteCarloBad_50[vect2, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="brown2")
mb5 = ggplot(df_monteCarloBad_100[vect2,], aes(y = df_monteCarloBad_100[vect2, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="orange")
mb6 = ggplot(df_monteCarloBad_500[vect2,], aes(y = df_monteCarloBad_500[vect2, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu2") + geom_boxplot(fill="aquamarine")
mb7 = ggplot(df_monteCarloBad_50[vect3,], aes(y = df_monteCarloBad_50[vect3, 3])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="brown2")
mb8 = ggplot(df_monteCarloBad_100[vect3,], aes(y = df_monteCarloBad_100[vect3, 3])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="orange")
mb9 = ggplot(df_monteCarloBad_500[vect3,], aes(y = df_monteCarloBad_500[vect3, 3])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour mu3") + geom_boxplot(fill="aquamarine")
ggarrange(mb1, mb2, mb3, mb4, mb5, mb6, mb7, mb8, mb9, ncol = 3, nrow = 3)


# Pour les \sigma

sb1 = ggplot(df_monteCarloBad_50[vect1,], aes(y = df_monteCarloBad_50[vect1, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="brown2")
sb2 = ggplot(df_monteCarloBad_100[vect1,], aes(y = df_monteCarloBad_100[vect1, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="orange")
sb3 = ggplot(df_monteCarloBad_500[vect1,], aes(y = df_monteCarloBad_500[vect1, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma1") + geom_boxplot(fill="aquamarine")
sb4 = ggplot(df_monteCarloBad_50[vect2,], aes(y = df_monteCarloBad_50[vect2, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="brown2")
sb5 = ggplot(df_monteCarloBad_100[vect2,], aes(y = df_monteCarloBad_100[vect2, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="orange")
sb6 = ggplot(df_monteCarloBad_500[vect2,], aes(y = df_monteCarloBad_500[vect2, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma2") + geom_boxplot(fill="aquamarine")
sb7 = ggplot(df_monteCarloBad_50[vect3,], aes(y = df_monteCarloBad_50[vect3, 4])) +
  xlab("n = 50") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="brown2")
sb8 = ggplot(df_monteCarloBad_100[vect3,], aes(y = df_monteCarloBad_100[vect3, 4])) +
  xlab("n = 100") + ylab("erreur absolue") +
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="orange")
sb9 = ggplot(df_monteCarloBad_500[vect3,], aes(y = df_monteCarloBad_500[vect3, 4])) +
  xlab("n = 500") + ylab("erreur absolue") + 
  ggtitle("Boxplot pour sigma3") + geom_boxplot(fill="aquamarine")
ggarrange(sb1, sb2, sb3, sb4, sb5, sb6, sb7, sb8, sb9, ncol = 3, nrow = 3)

# teste des fonctions simulation et algo_EM pour ce mélange de 4 gaussiennes
X4 = simulation(df_th4, 100)
print(df_th4)
algo_EM(data_init4, X4, 30)



# Chapitre III, étude d'un mélange de deux lois




bird_names2 = c("European Goldfinch", "Ring Ouzel")
mean_volume = c(38.0, 298.6)
sd_volume = c(9.1, 125.1)
nest_data2 = data.frame(bird_names2, mean_volume, sd_volume)
df <- random_species(nest_data2, 2)
df[,2] <- c(0.2878713, 0.7121287 )
df
data <- simulation(df, 500)
data

ggplot(data.frame(data), aes(x=data)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 100) +
  geom_density(alpha=.3, color = "purple", fill="darkcyan", size=1.5) +
  xlab('X') +
  ylab('Densité') +
  geom_vline(aes(xintercept = 40), col = "darkgreen", size = 1)+
  geom_vline(aes(xintercept = 300), col = "darkgreen", size = 1)+
  geom_vline(aes(xintercept = 10), col = "firebrick") + 
  geom_vline(aes(xintercept = 70), col = "firebrick") +
  geom_vline(aes(xintercept = 430), col = "firebrick") + 
  geom_vline(aes(xintercept = 170), col = "firebrick")


data_save <- data


param_init = data.frame(bird_names = c("European Goldfinch", "Ring Ouzel"),
                            alpha_init = c(0.5, 0.5), mean_init = c(40, 300),
                            sd_init = c(30, 130))


algo_EM(param_init, data, 10)




##### BOXPLOT #####

df_monteCarlo_100 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 100, 10)
df_monteCarlo_250 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 250, 10)
df_monteCarlo_500 = Monte_Carlo(good_dfTestTh3, good_dfTestInit3, 100, 500, 10)

# Pour les \alpha
#alpha1
df_Good1 = rbind(df_monteCarlo_100[vect1,],
                 df_monteCarlo_250[vect1,],
                 df_monteCarlo_500[vect1,])

df_Good1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Good1)[1] = "taille_echantillon"

a1 = qplot(df_Good1[,1], 
           df_Good1[,2], 
           data = df_Good1, 
           geom= "boxplot",
           fill = I("brown2")) + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[1])))
  
  

#alpha2
df_Good2 = rbind(df_monteCarlo_100[vect2,], 
                 df_monteCarlo_250[vect2,],
                 df_monteCarlo_500[vect2,])

df_Good2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Good2)[1] = "taille_echantillon"


a2 = qplot(df_Good2[,1],
           df_Good2[,2],
           data = df_Good2,
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") +
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[2])))

#alpha3
df_Good3 = rbind(df_monteCarlo_100[vect3,],
                 df_monteCarlo_250[vect3,],
                 df_monteCarlo_500[vect3,])

df_Good3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Good3)[1] = "taille_echantillon"


a3 = qplot(df_Good3[,1],
           df_Good3[,2],
           data = df_Good3,
           fill = I("aquamarine"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[3])))

# save with width = 800 and height = 1200


# Pour les \mu
# pour mu1
df_Good1 = rbind(df_monteCarlo_100[vect1,],
                 df_monteCarlo_250[vect1,],
                 df_monteCarlo_500[vect1,])

df_Good1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Good1)[1] = "taille_echantillon"


m1 = qplot(df_Good1[,1],
           df_Good1[,3],
           data = df_Good1, 
           geom= "boxplot",
           fill = I("brown2")) + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[1])))

#mu2
df_Good2 = rbind(df_monteCarlo_100[vect2,], 
                 df_monteCarlo_250[vect2,],
                 df_monteCarlo_500[vect2,])

df_Good2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Good2)[1] = "taille_echantillon"


m2 = qplot(df_Good2[,1],
           df_Good2[,3],
           data = df_Good2,
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[2])))

#mu3
df_Good3 = rbind(df_monteCarlo_100[vect3,],
                 df_monteCarlo_250[vect3,],
                 df_monteCarlo_500[vect3,])

df_Good3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Good3)[1] = "taille_echantillon"

m3 = qplot(df_Good3[,1], 
           df_Good3[,3], 
           data = df_Good3,
           fill = I("aquamarine"),
           geom= "boxplot")+
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[3])))

# save with width = 800 and height = 1200


# Pour les \sigma
# pour sigma1
df_Good1 = rbind(df_monteCarlo_100[vect1,],
                 df_monteCarlo_250[vect1,],
                 df_monteCarlo_500[vect1,])

df_Good1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Good)[1] = "taille_echantillon"

s1 = qplot(df_Good1[,1], 
           df_Good1[,4],
           data = df_Good1, 
           geom= "boxplot", 
           fill = I("brown2")) + 
  xlab("n") + ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[1])))

# Pour sigma2
df_Good2 = rbind(df_monteCarlo_100[vect2,],
                 df_monteCarlo_250[vect2,],
                 df_monteCarlo_500[vect2,])

df_Good2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Good2)[1] = "taille_echantillon"


s2 = qplot(df_Good2[,1],
           df_Good2[,4], 
           data = df_Good2, 
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[2])))

# Pour sigma3
df_Good3 = rbind(df_monteCarlo_100[vect3,], 
                 df_monteCarlo_250[vect3,],
                 df_monteCarlo_500[vect3,])

df_Good3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                 rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                 rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Good3)[1] = "taille_echantillon"

s3 = qplot(df_Good3[,1], 
           df_Good3[,4], 
           data = df_Good3, 
           fill = I("aquamarine"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[3])))


ggarrange(a1,a2,a3, ncol = 1, nrow = 3)
ggarrange(m1,m2,m3, ncol = 1, nrow = 3)
ggarrange(s1,s2,s3, ncol = 1, nrow = 3)

########################################################################
# Pour les \alpha
#alpha1
df_Bad1 = rbind(df_monteCarlo_100[vect1,],
                df_monteCarlo_250[vect1,],
                df_monteCarlo_500[vect1,])

df_Bad1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Bad1)[1] = "taille_echantillon"

a1 = qplot(df_Bad1[,1], 
           df_Bad1[,2], 
           data = df_Bad1, 
           geom= "boxplot",
           fill = I("brown2")) + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[1])))

#alpha2
df_Bad2 = rbind(df_monteCarlo_100[vect2,], 
                df_monteCarlo_250[vect2,],
                df_monteCarlo_500[vect2,])

df_Bad2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Bad2)[1] = "taille_echantillon"


a2 = qplot(df_Bad2[,1],
           df_Bad2[,2],
           data = df_Bad2,
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") +
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[2])))

#alpha3
df_Bad3 = rbind(df_monteCarlo_100[vect3,],
                df_monteCarlo_250[vect3,],
                df_monteCarlo_500[vect3,])

df_Bad3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Bad3)[1] = "taille_echantillon"


a3 = qplot(df_Bad3[,1],
           df_Bad3[,2],
           data = df_Bad3,
           fill = I("aquamarine"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", alpha[3])))

# save with width = 800 and height = 1200


# Pour les \mu
# pour mu1
df_Bad1 = rbind(df_monteCarlo_100[vect1,],
                df_monteCarlo_250[vect1,],
                df_monteCarlo_500[vect1,])

df_Bad1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Bad1)[1] = "taille_echantillon"


m1 = qplot(df_Bad1[,1],
           df_Bad1[,3],
           data = df_Bad1, 
           geom= "boxplot",
           fill = I("brown2")) + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[1])))

#mu2
df_Bad2 = rbind(df_monteCarlo_100[vect2,], 
                df_monteCarlo_250[vect2,],
                df_monteCarlo_500[vect2,])

df_Bad2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Bad2)[1] = "taille_echantillon"


m2 = qplot(df_Bad2[,1],
           df_Bad2[,3],
           data = df_Bad2,
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[2])))

#mu3
df_Bad3 = rbind(df_monteCarlo_100[vect3,],
                df_monteCarlo_250[vect3,],
                df_monteCarlo_500[vect3,])

df_Bad3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Bad3)[1] = "taille_echantillon"

m3 = qplot(df_Bad3[,1], 
           df_Bad3[,3], 
           data = df_Bad3,
           fill = I("aquamarine"),
           geom= "boxplot")+
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", mu[3])))

# save with width = 800 and height = 1200


# Pour les \sigma
# pour sigma1
df_Bad1 = rbind(df_monteCarlo_100[vect1,],
                df_monteCarlo_250[vect1,],
                df_monteCarlo_500[vect1,])

df_Bad1[,1] = c(rep("100", dim(df_monteCarlo_100[vect1,])[1]),
                rep("250", dim(df_monteCarlo_250[vect1,])[1]),
                rep("500", dim(df_monteCarlo_500[vect1,])[1]))

colnames(df_Bad1)[1] = "taille_echantillon"

s1 = qplot(df_Bad1[,1], 
           df_Bad1[,4],
           data = df_Bad1, 
           geom= "boxplot", 
           fill = I("brown2")) + 
  xlab("n") + ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[1])))

# Pour sigma2
df_Bad2 = rbind(df_monteCarlo_100[vect2,],
                df_monteCarlo_250[vect2,],
                df_monteCarlo_500[vect2,])

df_Bad2[,1] = c(rep("100", dim(df_monteCarlo_100[vect2,])[1]),
                rep("250", dim(df_monteCarlo_250[vect2,])[1]),
                rep("500", dim(df_monteCarlo_500[vect2,])[1]))

colnames(df_Bad2)[1] = "taille_echantillon"


s2 = qplot(df_Bad2[,1],
           df_Bad2[,4], 
           data = df_Bad2, 
           fill = I("orange"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[2])))

# Pour sigma3
df_Bad3 = rbind(df_monteCarlo_100[vect3,], 
                df_monteCarlo_250[vect3,],
                df_monteCarlo_500[vect3,])

df_Bad3[,1] = c(rep("100", dim(df_monteCarlo_100[vect3,])[1]),
                rep("250", dim(df_monteCarlo_250[vect3,])[1]),
                rep("500", dim(df_monteCarlo_500[vect3,])[1]))

colnames(df_Bad3)[1] = "taille_echantillon"

s3 = qplot(df_Bad3[,1], 
           df_Bad3[,4], 
           data = df_Bad3, 
           fill = I("aquamarine"),
           geom= "boxplot") + 
  xlab("n") + 
  ylab("erreur absolue") +
  ggtitle(expression(paste("Boxplot des erreurs pour ", sigma[3])))


ggarrange(a1,a2,a3, ncol = 1, nrow = 3)
ggarrange(m1,m2,m3, ncol = 1, nrow = 3)
ggarrange(s1,s2,s3, ncol = 1, nrow = 3)
