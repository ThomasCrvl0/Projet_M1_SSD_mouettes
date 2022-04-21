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

# Fonction retournant la liste composée de n couples (moyenne, variance)
# tirés aléatoirement parmi les 13 espèces d'oiseaux
# Ces n couples serviront à tirer aléatoirement 100 valeurs 
# issues d'un mélange gaussien 
nb_species = function(n){
  index_species = sample(1:dim(nest_data)[1], n, replace = FALSE)
  listOfparam = list()
  for (i in index_species){
    l_param = list(as.vector(t(nest_data[i,c(2,3)])))
    listOfparam = append(listOfparam, l_param)
  }
  return(listOfparam)
}


GetRealParam = function(n){
  if (n == 2){
    real_mu1 = listOfparam[[1]][1]
    real_sd1 = listOfparam[[1]][2]
    real_mu2 = listOfparam[[2]][1]
    real_sd2 = listOfparam[[2]][2]
    real_param = data.frame(real_mean = c(real_mu1, real_mu2),
                            real_sd = c(real_sd1, real_sd2))
  } else if (n == 3){
    real_mu1 = listOfparam[[1]][1]
    real_sd1 = listOfparam[[1]][2]
    real_mu2 = listOfparam[[2]][1]
    real_sd2 = listOfparam[[2]][2]
    real_mu3 = listOfparam[[3]][1]
    real_sd3 = listOfparam[[3]][2]
    real_param = data.frame(real_mean = c(real_mu1, real_mu2, real_mu3),
                            real_sd = c(real_sd1, real_sd2, real_sd3))
  } else {
    real_mu1 = listOfparam[[1]][1]
    real_sd1 = listOfparam[[1]][2]
    real_mu2 = listOfparam[[2]][1]
    real_sd2 = listOfparam[[2]][2]
    real_mu3 = listOfparam[[3]][1]
    real_sd3 = listOfparam[[3]][2]
    real_mu4 = listOfparam[[3]][1]
    real_sd4 = listOfparam[[3]][2]
    real_param = data.frame(real_mean = c(real_mu1, real_mu2, real_mu3, real_mu4),
                            real_sd = c(real_sd1, real_sd2, real_sd3, real_sd4))
  }
  return(real_param)
}


matrixOfMixture = function(n){
  X = NULL
  X = matrix(nrow=100,ncol=1)
  ### Boucle permettant de tirer 100 valeurs issues
  ### d’un mélange gaussien :
  if (n == 2){
    for (i in 1:100) {
      Z = rbinom(1,1,0.5) # choix de la loi par tirage de Bernoulli
      if (Z == 1) {
        mu1 = listOfparam[[1]][1]
        sd1 = listOfparam[[1]][2]
        X[i] = rnorm(1,mu1,sd1)
      }else {
        mu2 = listOfparam[[2]][1]
        sd2 = listOfparam[[2]][2]
        X[i] = rnorm(1,mu2,sd2)
      }
    }
  }else if(n == 3){
    for (i in 1:100) {
      Z = sample(1:3, 1)
      if (Z == 1) {
        mu1 = listOfparam[[1]][1]
        sd1 = listOfparam[[1]][2]
        X[i] = rnorm(1,mu1,sd1)
      }else if(Z == 2){
        mu2 = listOfparam[[2]][1]
        sd2 = listOfparam[[2]][2]
        X[i] = rnorm(1,mu2,sd2)
      }else{
        mu3 = listOfparam[[3]][1]
        sd3 = listOfparam[[3]][2]
        X[i] = rnorm(1,mu3,sd3)
      }
    }
  }else{
    for (i in 1:100) {
      Z = sample(1:4, 1)
      if (Z == 1) {
        mu1 = listOfparam[[1]][1]
        sd1 = listOfparam[[1]][2]
        X[i] = rnorm(1,mu1,sd1)
      }else if(Z == 2){
        mu2 = listOfparam[[2]][1]
        sd2 = listOfparam[[2]][2]
        X[i] = rnorm(1,mu2,sd2)
      }else if(Z == 3){
        mu3 = listOfparam[[3]][1]
        sd3 = listOfparam[[3]][2]
        X[i] = rnorm(1,mu3,sd3)
      }else{
        mu4 = listOfparam[[4]][1]
        sd4 = listOfparam[[4]][2]
        X[i] = rnorm(1,mu4,sd4)
      }
    }
  }
  
  return(X)
}


# on fera K itérations de l'algorithme
simulation = function(n,K){
  df_th_param = GetRealParam(n)
  if (n == 2){
    cste = sample(c(-15,15), 4, replace = TRUE)
    ## Définition arbitraire des valeurs initiales des paramètres :
    lambda1 = 0.2
    lambda2 = 0.8
    mu1 = df_th_param$real_mean[1] + cste[1]
    mu2 = df_th_param$real_mean[2] + cste[2]
    sigma1 = abs(df_th_param$real_sd[1] + cste[3])
    sigma2 = abs(df_th_param$real_sd[2] + cste[4])
    for (i in 1:K) {
      ## Application de la formule (28) :
      vrais1 = lambda1*dnorm(X,mean=mu1,sd=sigma1)
      vrais2 = lambda2*dnorm(X,mean=mu2,sd=sigma2)
      vrais12 = vrais1/(vrais1 + vrais2) # probas a posteriori p_{i,1}
      vrais22 = vrais2/(vrais1 + vrais2) # probas a posteriori p_{i,2}
      ## Mise à jour de lambda1 et lambda2 = P(Z=1 | X,Thteta) et P(Z=2 | X,Thteta):
      lambda1 = mean(vrais12)
      lambda2 = 1 - lambda1
      ## Mise à jour de mu1 et mu2 :
      mu1 = sum(vrais12*X)/sum(vrais12)
      mu2 = sum(vrais22*X)/sum(vrais22)
      ## Mise à jour de sigma1 et sigma2 :
      sigma1 = sqrt(sum(vrais12*(X-mu1)^2)/(sum(vrais12)))
      sigma2 = sqrt(sum(vrais22*(X-mu2)^2)/(sum(vrais22)))
    }
    new_df = data.frame(estimate_mean = c(mu1, mu2),
                        estimate_sd = c(sigma1, sigma2))
  }else if(n == 3){
    cste = sample(c(-30,30), 6, replace = TRUE)
    ## Définition arbitraire des valeurs initiales des paramètres :
    lambda1 = 0.2
    lambda2 = 0.5
    lambda3 = 0.3
    mu1 = df_th_param$real_mean[1] + cste[1]
    mu2 = df_th_param$real_mean[2] + cste[2]
    mu3 = df_th_param$real_mean[3] + cste[3]
    sigma1 = abs(df_th_param$real_sd[1] + cste[4])
    sigma2 = abs(df_th_param$real_sd[2] + cste[5])
    sigma3 = abs(df_th_param$real_sd[3] + cste[6])
    for (i in 1:K) {
      ## Application de la formule (28) :
      vrais1 = lambda1*dnorm(X,mean=mu1,sd=sigma1)
      vrais2 = lambda2*dnorm(X,mean=mu2,sd=sigma2)
      vrais3 = lambda3*dnorm(X,mean=mu3,sd=sigma3)
      vrais13 = vrais1/(vrais1 + vrais2 + vrais3) # probas a posteriori p_{i,1}
      vrais23 = vrais2/(vrais1 + vrais2 + vrais3) # probas a posteriori p_{i,2}
      vrais33 = vrais3/(vrais1 + vrais2 + vrais3) # probas a posteriori p_{i,3}
      ## Mise à jour de lambda_j pour j = {1,2,3}:
      lambda1 = mean(vrais13)
      lambda2 = mean(vrais23)
      lambda3 = 1 - lambda1 - lambda2
      ## Mise à jour de mu1, mu2 et mu3:
      mu1 = sum(vrais13*X)/sum(vrais13)
      mu2 = sum(vrais23*X)/sum(vrais23)
      mu3 = sum(vrais33*X)/sum(vrais33)
      ## Mise à jour de sigma1, sigma2 et sigma3:
      sigma1 = sqrt(sum(vrais13*(X-mu1)^2)/(sum(vrais13)))
      sigma2 = sqrt(sum(vrais23*(X-mu2)^2)/(sum(vrais23)))
      sigma3 = sqrt(sum(vrais33*(X-mu3)^2)/(sum(vrais33)))
    }
    new_df = data.frame(estimate_mean = c(mu1, mu2, mu3),
                        estimate_sd = c(sigma1, sigma2, sigma3))
  }else{
    cste = sample(c(-30,30), 8, replace = TRUE)
    ## Définition arbitraire des valeurs initiales des paramètres :
    lambda1 = 0.2
    lambda2 = 0.4
    lambda3 = 0.3
    lambda4 = 0.1
    mu1 = df_th_param$real_mean[1] + cste[1]
    mu2 = df_th_param$real_mean[2] + cste[2]
    mu3 = df_th_param$real_mean[3] + cste[3]
    mu4 = df_th_param$real_mean[4] + cste[4]
    sigma1 = abs(df_th_param$real_sd[1] + cste[5])
    sigma2 = abs(df_th_param$real_sd[2] + cste[6])
    sigma3 = abs(df_th_param$real_sd[3] + cste[7])
    sigma4 = abs(df_th_param$real_sd[4] + cste[8])
    for (i in 1:K) {
      ## Application de la formule (28) :
      vrais1 = lambda1*dnorm(X,mean=mu1,sd=sigma1)
      vrais2 = lambda2*dnorm(X,mean=mu2,sd=sigma2)
      vrais3 = lambda3*dnorm(X,mean=mu3,sd=sigma3)
      vrais4 = lambda4*dnorm(X,mean=mu4,sd=sigma4)
      # probas a posteriori p_{i,1}
      vrais14 = vrais1/(vrais1 + vrais2 + vrais3 + vrais4)
      # probas a posteriori p_{i,2}
      vrais24 = vrais2/(vrais1 + vrais2 + vrais3 + vrais4)
      # probas a posteriori p_{i,3}
      vrais34 = vrais3/(vrais1 + vrais2 + vrais3 + vrais4)
      # probas a posteriori p_{i,4}
      vrais44 = vrais4/(vrais1 + vrais2 + vrais3 + vrais4)
      ## Mise à jour de lambda_j pour j = {1,2,3,4}:
      lambda1 = mean(vrais14)
      lambda2 = mean(vrais24)
      lambda3 = mean(vrais34)
      lambda4 = 1 - lambda1 - lambda2 - lambda3
      ## Mise à jour de mu1, mu2, mu3 et mu4:
      mu1 = sum(vrais14*X)/sum(vrais14)
      mu2 = sum(vrais24*X)/sum(vrais24)
      mu3 = sum(vrais34*X)/sum(vrais34)
      mu4 = sum(vrais44*X)/sum(vrais44)
      ## Mise à jour de sigma1, sigma2, sigma3, sigma4 :
      sigma1 = sqrt(sum(vrais14*(X-mu1)^2)/(sum(vrais14)))
      sigma2 = sqrt(sum(vrais24*(X-mu2)^2)/(sum(vrais24)))
      sigma3 = sqrt(sum(vrais34*(X-mu3)^2)/(sum(vrais34)))
      sigma4 = sqrt(sum(vrais44*(X-mu4)^2)/(sum(vrais44)))
    }
    new_df = data.frame(estimate_mean = c(mu1, mu2, mu3, mu4),
                        estimate_sd = c(sigma1, sigma2, sigma3, sigma4))
  }
  return(new_df)
}


# Execution du code
n = 2; K = 30
listOfparam = nb_species(n)
X = matrixOfMixture(n)
print(GetRealParam(n))
simulation(n, K)