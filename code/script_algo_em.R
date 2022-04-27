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
df_th = GetRealParam(c("European Goldfinch", "Ring Ouzel"))

# Dataframe des données d'initialisations choisies par l'utilisateur
data_init = data.frame(species_chosen = c("European Goldfinch", "Ring Ouzel")
                       ,alpha_init = c(0.2,0.8), mean_init = c(50, 280),
                       sd_init = c(11, 130))

simulation = function(data, n=100){
  X = NULL
  X = matrix(nrow=n,ncol=1)
  J = dim(data)[1]
  for(i in 1:n){
    Z = runif(1)
    if (Z <= data$alpha_init[1]){
      X[i] = rnorm(1, data$mean_init[1], data$sd_init[1])
    }else{
      k = 1
      l = 2
      Bool = FALSE
      vec_alpha = data$alpha_init
      cumul_alpha = cumsum(vec_alpha)
      while((Bool == FALSE) & (k < J) & (l < J+1)){
        if((cumul_alpha[k]<=Z) & (cumul_alpha[l]>=Z)){
          Bool = TRUE
          param_index = l
        }
        k = k+1
        l = l+1
      }
      X[i] = rnorm(1, data$mean_init[param_index], data$sd_init[param_index])
    }
  }
  return(X)
}
