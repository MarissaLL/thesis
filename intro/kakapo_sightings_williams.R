# This script used for data entry and correction. The final output is a datafile which is imported into plot_kakapo_map.R

# Historical sites from Williams et al. 1956
introduced_NI <- tibble(number = c(1, 19),
                     site = c("Little Barrier Island", "Kapiti Island"),
                     latitude = c(-36.20252,-40.85372),
                     longitude = c(175.08098, 174.91641))

subfossil_NI <- tibble(number = c(3:6, 13:18),
                       site = c("Frankton", "Matiri?", "Karamu","Waitomo", "Coonoor", "Akitio", "Pahiatua", "Waewaepa", "Mauriceville", "Martinborough"),
                       latitude = c(-37.78983, -37.80406, -37.88580, -38.26135, -40.43477, -40.60408, -40.45444, -40.36824, -40.77592,-41.21853),
                       longitude = c(175.25086, 175.39401,  175.13785,175.10354, 176.11315, 176.40962, 175.84105, 176.14847, 175.69988, 175.45948))



hs_williams_NI <- tibble(number = c(2, 7:12),
                                       site = c("Hunua Range", "Upper Wanganui River", "Raurimu - Wanganui River", "Kaimanawa Ranges", "Kaimanawa Ranges", "Ngaruroro River", 
                                                "Huiarau Ranges"), 
                                       latitude = c(-37.07506, -39.31618, -39.19943, -39.046919, -39.24182, -39.436193, -38.51659), 
                                       longitude = c(175.18006, 174.95317, 175.20118, 176.121826, 175.88011, 176.301126, 177.24080)) 




subfossil_SI <- tibble(number = c(22:32, 37),
                       site = c("Lake Grassmere",  "Pyramid Valley", "Broken River",  "Sumner", "Mt Somers", "Pareora", "Timaru", "Castle Rock", "Owaka", "Forest Hills?", "Takaka", "Karamea"),
                       latitude = c(-41.737984, -42.96979, -43.20068, -43.56929, -43.67270, -44.42204, -44.40698, -45.80483, -46.45294, -45.963561, -40.85085, -41.24817),
                       longitude = c(174.196301, 172.60020, 171.78815, 172.75959, 171.36025, 171.03515, 171.24926, 168.29354, 169.65714, 168.249435, 172.80703, 172.11219))





hs_williams_SI_Nelson_Nwestland <-  tibble(number = c(20:21, 33:36, 38:47),
                                               site = c( "Mt Stokes", "Maori Bay", "Gouland Downs", "Gouland Downs", "Leslie Valley and Kakapo Valley",  "Mt Arthur",  "Mt Owen", 
                                                        "Westport. area", "Moonlight", "Grey-Inangahua Saddle", "Robinson River", "Bold Head", "Ross", "Rangitata", "Okarito",  "Bruce Bay"),
                                               latitude = c(-41.09042, -41.17182, -40.896257, -40.927132, -41.25500, -41.21791, -41.55200, -41.84160, -42.27030,  -42.30900, -42.49936, -42.94051,-42.92452, -43.57278, -43.22329, -43.60730),
                                               longitude = c(174.10196, 173.83205, 172.302876 , 172.366047, 172.41847, 172.68152, 172.54131, 171.76660, 171.45070, 171.98408, 172.08411, 170.70447, 170.83089, 170.76292, 
                                                             170.16169, 169.59137 ))



hs_williams_SI_Bruce_Wanaka <- tibble(number = c(48:58), 
                                                    site = c("Paringa", "Okuru-Haast divide", "Haast-Clarke junction", "Haast Pass", "Burke Flat", "Jackson Bay", "Jackson Bay", "Jackson Bay", "Big Bay",  
                                                             "Five Fingers Range", "Red Pyke River"),
                                                    latitude = c(-43.80615, -44.09447, -43.94690, -44.10721,  -44.01973, -44.02068, -44.071554, -44.076734, -40.85826, -44.41313, -44.31252),
                                                    longitude = c(169.51717, 169.24160, 169.44979, 169.35481, 169.36837, 168.64544, 168.787079, 168.634644, 172.14202, 168.43120, 168.31141)) 

introduced_SI <- tibble(number = c(100:104, 111),
                        site = c("Long Island", "Resolution Island", "Parrot Island", "Anchor Island",  "Indian Island", "Christmas Village"),
                        latitude = c(-45.76146, -45.66042, -45.70720, -45.75877, -45.77818, -46.75111),
                        longitude = c(166.69701, 166.63798, 166.53330, 166.51186, 166.58723, 167.98087))

                                                                                                    
hs_williams_SI_Martin_Milford <- tibble(number = c(59:69), 
                                                      site = c("Martin Bay and Sara Hills", "Martin Bay area", "Lake McKerrow", "Hollyford Valley", "Tutoko Valley", "Mitre Peak", "Lake Ada", "Milford Sound", "Cleddau Valley", 
                                                               "Homer Tunnel", "Upper Hollyford Valley"),
                                                      latitude = c(-44.36544, -44.40400, -44.49420, -44.59310, -44.61335,-44.63924, -44.70853, -44.66523, -44.71478, -44.76419, -44.80224),
                                                      longitude = c(168.05318, 167.96724, 168.07336, 168.12740, 167.98148, 167.83063, 167.85760, 167.92564, 167.95699, 167.98124,  168.02592))


hs_williams_SI_Milford_TeAnau <- tibble(number = c(70:91), 
                                                      site = c("Franklin Range", "Wild Native River", "Lake Beddoes", "Glaisnock-Wild Native divide", "Glaisnock Valley and Nitz Creek", "Henry Saddle", "Stillwater River", "George Sound", 
                                                               "Vigil Peak and Ethne Saddle area", "East of Whitewater River", "Leslie Clearing", "Caswell Sound", "Mt Pisgah", "Birch Creek", "Turret Peaks", "Mid Burn", 
                                                               "Middle Fiord area", "Ettrick Burn area", "Takahe Valley", "Coronation Peak area", "Wilmot Pass area", "Wilmot Pass area"),
                                                      latitude = c(-44.87648, -44.88544, -44.90663, -44.93638,-44.94988, -45.00775, -45.02475, -44.99472,  -45.01163, -44.94570, -45.01696, -45.04799, 
                                                                   -45.10971, -45.08275, -45.12206, -45.14177, -45.17480, -45.25214, -45.28912, -45.22050, -45.50817,-45.54467),
                                                      longitude = c(167.65616, 167.56891, 167.54603, 167.58919, 167.71837, 167.49724, 167.41679, 167.41385, 167.39694, 167.32872,  167.30862, 167.30884, 
                                                                    167.49696, 167.71189, 167.75549, 167.74318, 167.66103, 167.62674,  167.65743, 167.26616, 167.19254, 167.19124))


hs_williams_SI_Doubtful_Hauroko <- tibble(number = c(92:99, 105:110),
                                                        site = c("Vancouver Arm - Hall Arm divide", "Mackenzie Flat area", "Roa Stream", "Roa Stream", "Green Point?", "Cooper Island", "Wet Jacket Arm",  "Breaksea Sound area",   
                                                                 "Cascade Cove", "Pickersgill Harbour", "Tableland", "Northport", "Chalky Sound area", "Preservation Inlet area"),
                                                        latitude = c(-45.47375, -45.57266, -45.70489,-45.71620, -45.72478, -45.73272,  -45.65171, -45.58618,-45.816986,-45.797481, -45.930617,-45.971674, -46.029866, -46.090019),
                                                        longitude = c(167.03339, 167.13847, 167.02051, 167.04770, 166.92970, 166.82844, 166.80881, 166.81976, 166.565866, 166.572647, 166.559258, 166.576767, 166.615219,  166.721821))

hs_williams_stewartI <- tibble(number = c(112:113),
                                             site = c( "Port Pegasus", "Port Pegasus-Seal Point"),
                                             latitude = c(-47.158556, -47.176995),
                                             longitude = c(167.723808, 167.822685))


historical_sites_williams <- bind_rows(hs_williams_NI,  hs_williams_SI_Nelson_Nwestland, hs_williams_SI_Bruce_Wanaka, hs_williams_SI_Martin_Milford,
                                       hs_williams_SI_Milford_TeAnau, hs_williams_SI_Doubtful_Hauroko, hs_williams_stewartI) %>% 
  mutate(type = "Historical sighting") %>% 
  rename(unrounded_lat = latitude, unrounded_long = longitude) %>% 
  mutate(latitude = round(unrounded_lat, digits = 2),
         longitude = round(unrounded_long, digits = 2))

historical_introductions <- bind_rows(introduced_NI, introduced_SI) %>% 
  mutate(type = "Introduction") %>% 
  rename(unrounded_lat = latitude, unrounded_long = longitude) %>% 
  mutate(latitude = round(unrounded_lat, digits = 2),
         longitude = round(unrounded_long, digits = 2))

subfossil <- bind_rows(subfossil_NI, subfossil_SI) %>% 
  mutate(type = "Subfossil remains")%>% 
  rename(unrounded_lat = latitude, unrounded_long = longitude) %>% 
  mutate(latitude = round(unrounded_lat, digits = 2),
         longitude = round(unrounded_long, digits = 2))


all <- bind_rows(historical_sites_williams, historical_introductions, subfossil) %>% 
  arrange(number) %>% 
  select(-c(unrounded_lat, unrounded_long))

# write_delim(historical_sites_williams, "~/projects/kakapo-genomics/data/historic_kakapo_sightings_williams56.csv", delim = ",")
# write_delim(all, "~/projects/kakapo-genomics/data/all_historic_kakapo_williams56.csv", delim = ",")


# 
# # Test plot #
# # Load in map from maps package
# nz_map <- map_data("nz")
# 
# 
# # Kakapo colours from manu: dark green "#7D9D33" light green "#CED38C" yellow "#DCC949" grey-brown "#BCA888" red-brown"#CD8862" dark-brown "#775B24"
# # Plot the map
# ggplot() + 
#   geom_polygon(data = nz_map, 
#                mapping = aes(x=long, 
#                              y = lat, 
#                              group = group), 
#                fill = NA, colour = 'black' ) + 
#   coord_fixed(1.4) + 
#   geom_point(data =  all,
#              mapping = aes(x = longitude,
#                            y = latitude, colour = type),
#              size = 3.5, alpha = 0.7) +
#   scale_color_manual(values = c( "#CED38C",  "#CD8862", "#7D9D33")) + 
#   theme_classic() +
#   theme(axis.line.x.top = element_line(),
#         axis.line.y.right = element_line(),
#         panel.background = element_rect(fill = NA,
#                                         colour = "black",
#                                         size = 0.7 ),
#         axis.title = element_text(size = 18, margin = margin(r = 12)),
#         axis.text = element_text(size = 12),
#         legend.title = element_blank()) +
#   labs(x = 'Longitude',
#        y = 'Latitude')#+
#   # geom_text(data = done_so_far, 
#   #                     mapping = aes(x = longitude,
#   #                                   y = latitude,
#   #                                   label = number),
#   #                     hjust = -1,
#   #                     size = 5)
# 
# 
# 

















# OLD VER IN CASE I MESS UP THE NEW VER.
# 
# # Historical sites from Williams et al. 1956
# historical_sites_williams_NI <- tibble(number = c(1:19),
#                                        site = c("Little Barrier Island", "Hunua Range", "Frankton", "Matiri", "Karamu", "Waitomo", "Upper Wanganui River", "Raurimu - Wanganui River", "Kaimanawa Ranges", "Kaimanawa Ranges", "Ngaruroro River", 
#                                                 "Huiarau Ranges", "Coonoor", "Akitio", "Pahiatua", "Waewaepa", "Mauriceville", "Martinborough", "Kapiti Island"))
# historical_sites_williams_SI_E <- tibble(number = c(20:31),
#                                          site = c("Mt Stokes", "Maori Bay",  "Lake Grassmere",  "Pyramid Valley", "Broken River",  "Sumner", "Mt Somers", "Pareora", "Timaru", "Castle Rock", "Owaka", "Forest Hills"))
# 
# # Should only go to 47?                                    
# historical_sites_williams_SI_nelson <-  tibble(number = c(32:49),
#                                                site = c("Takaka", "Gouland Downs", "Gouland Downs", "Leslie Valley", "Kakapo Valley ND",  "Mt Arthur", "Karamea", "Mt Owen", 
#                                                         "Northern Westland", "Westport. area", "Moonlight", "Grey-Inangahua Saddle", "Robinson River", "Bold Head ND", "Ross", "Rangitata", "Okarito",  "Bruce Bay")) 
# historical_sites_williams_SI_Bruce_Wanaka <- tibble(number = c(48:58), 
#                                                     site = c("Paringa", "Okuru-Haast divide", "Haast-Clarke junction", "Haast Pass", "Burke Bay", "Jackson Bay", "Jackson Bay ND", "Jackson Bay", "Big Bay",  
#                                                              "Five Fingers Range", "Red Pike River"))
# 
# ## Should only go to 69?                                                                                                    
# historical_sites_williams_SI_Martin_Milford <- tibble(number = c(59:70), 
#                                                       site = c("Martin Bay", "Sara Hills", "Martin Bay area", "Lake McKerrow", "Hollyford Valley", "Tutoko Valley", "Mitre Peak", "Lake Ada", "Milford
#                                                                                                      Sound", "Cleddau Valley", "Homer Tunnel", "Upper Hollyford Valley"))
# historical_sites_williams_SI_Milford_TeAnau <- tibble(number = c(70:91), 
#                                                       site = c("Franklin Range", "Wild Native R.", "Lake Beddoes", "Glaisnock-Wild Native divide", "Glaisnock Valley and Nitz Creek", "Henry Saddle", "Stillwater River", "George Sound", 
#                                                                "Vigil Peak and Ethne Saddle area", "East of Whitewater R.", "Leslie Clearing", "Caswell Sound", "Mt Pisgah", "Birch Creek", "Turret Peaks", "Mid Burn", "Middle Fiord area", 
#                                                                "Ettrick Burn area", "Takahe Valley", "Coronation Peak area", "Wilmot Pass area", "Wilmot Pass area"))
# 
# 
# historical_sites_williams_SI_Doubtful_Hauroko <- tibble(number = c(92:110),
#                                                         site = c("Vancouver Arm - Hall Arm divide", "Mackenzie Bay area", "Roa Stream", "Roa Stream", "Green Point", "Cooper Island", "Wet Jacket Arm",  "Breaksea Sound area",  "Long Island", 
#                                                                  "Resolution Island", "Parrot Island", "Anchor Island",  "Indian Island", "Cascade Cove", "Pickersgill Harbour", "Tableland", "Northport", "Chalky Sound area", "Preservation Inlet area"))
# 
# historical_sites_williams_stewartI <- tibble(number = c(111:113),
#                                              site = c("Christmas Village", "Port Pegasus", "Port Pegasus-Seal Point"))
# 
# 
# historical_sites_williams <- bind_rows(historical_sites_williams_NI, historical_sites_williams_SI_E, historical_sites_williams_SI_nelson, historical_sites_williams_SI_Bruce_Wanaka, historical_sites_williams_SI_Martin_Milford,
#                                        historical_sites_williams_SI_Doubtful_Hauroko, historical_sites_williams_SI_Milford_TeAnau, historical_sites_williams_stewartI)
