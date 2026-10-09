###############################################################################
# Escala de Actitudes hacia el Narcotráfico (ATN) - Reynoso et al. (2018)
# Reconstrucción SINTÉTICA a partir de estadísticos publicados + re-análisis
#
# QUÉ ES: simulación "moment-matching" (bootstrap paramétrico) condicionada a
#         lo publicado. Prueba coherencia interna y reproducibilidad analítica.
#
# NUMERACIÓN: i01..i24 = numeración de la Tabla 1 (versión de 24 ítems).
# Equivalencia Anexo 1 (17 ítems) -> Tabla 1:
#   Anexo:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
#   Tabla1: 1  4  5  6  7  9 10 11 12 13 14 15 20 21 22 23 24
#
# SUPUESTOS:
#  A1. Las cargas de la Tabla 1 (PAF + Equamax, ortogonal) se usan como modelo
#      de correlaciones PEARSON entre ítems ordinales: R = L L' + diag(1-h2).
#  A2. Dentro de cada familia (Rechazo / Apoyo / Predisposición) los ítems
#      comparten marginal. Marginales de ítems eliminados: supuestas.
#  A3. Cópula gaussiana; ítem = redondeo+recorte a 1..5 de una latente normal.
#  A4. Sexo y nivel educativo independientes (la conjunta no se reporta).
###############################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))   

pkgs <- c("psych", "lavaan", "GPArotation")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
invisible(lapply(pkgs, library, character.only = TRUE))
set.seed(20181031)

## ---------------------------------------------------------------------------
## 0. DATOS PUBLICADOS
## N: Tamaño de muestra
## n_bach / n_lic: cantidad de personas con ese nivel de estudio
## p_hombre: probabilidad de ser hombre
## p_bach: probabilidad de ser de bachillerato
## ---------------------------------------------------------------------------
N <- 2356; n_bach <- 568; n_lic <- 1788
p_hombre <- 0.494; p_bach <- n_bach / N
nm <- function(v) sprintf("i%02d", v)

#Tabla 1 del paper: cargas factoriales + Equamax (24 ítems, 3 factores)
tab1 <- read.table(header = TRUE, text = "
item f1 f2 f3
20  .747 -.064 -.129
14  .682 -.020 -.185
15  .670 -.068 -.137
23  .658 -.049 -.111
22  .619 -.224  .000
5   .603 -.022 -.186
10  .600 -.035 -.112
7   .576 -.098 -.046
12  .513 -.050 -.251
3   .378 -.054 -.216
13 -.085  .647  .171
6  -.153  .543  .139
11 -.176  .528  .325
9  -.005  .485  .194
21 -.329  .455  .187
17 -.143  .438  .366
19 -.003  .424 -.103
8   .021  .365  .059
2   .039  .195  .082
4  -.033  .315  .572
1  -.050  .354  .546
16  .479  .044 -.480
18  .367  .076 -.457
24 -.124  .327  .452
")
#Ordenar items en orden de aparición y no por carga
tab1 <- tab1[order(tab1$item), ]
#Convierte tabla en la matriz Lambda
Lam  <- as.matrix(tab1[, c("f1", "f2", "f3")])
rownames(Lam) <- nm(tab1$item)

#Escala final. (El texto del paper lista "17" en el Factor 2);
# la Tabla 1 y el Anexo indican que es el ítem 9.
F1 <- c(5, 7, 10, 12, 14, 15, 20, 22, 23)   #Rechazo (9)
F2 <- c(6, 9, 11, 13, 21)                   #Apoyo (5)
F3 <- c(1, 4, 24)                           #Predisposición (3)
ITEMS17 <- c(F1, F2, F3)
ITEMS18 <- c(F1, F2, 19, F3)                #modelos 1 y 2 del AFC

#Familia de marginal de cada ítem: 1 = Rechazo, 2 = Apoyo, 3 = Predisp.
fam <- rep(NA_integer_, 24)
fam[c(3, 5, 7, 10, 12, 14, 15, 16, 18, 20, 22, 23)]  <- 1
fam[c(2, 6, 8, 9, 11, 13, 17, 19, 21)]               <- 2
fam[c(1, 4, 24)]                                     <- 3

#Tablas 3, 5, 6. Orden: general, hombres, mujeres, bachillerato, licenciatura
tg <- list(
  #Factor 1
  list(mu = c(3.65, 3.62, 3.70, 3.25, 3.79), sd = c(.945, .978, .912, .973, .900)),
  #Factor 2
  list(mu = c(1.35, 1.48, 1.24, 1.66, 1.26), sd = c(.656, .749, .526, .840, .554)),
  #Factor 3
  list(mu = c(1.70, 1.79, 1.63, 1.95, 1.63), sd = c(.723, .764, .673, .819, .673)))
tot_pub <- c(gen = 1.98, H = 2.05, M = 1.92, B = 2.32, L = 1.87)

## ---------------------------------------------------------------------------
## 1. AUDITORÍA DE COHERENCIA INTERNA:
##    revisar si los números publicados por los autores tienen sentido teórico
## ---------------------------------------------------------------------------
cat("\n==== 1. AUDITORÍA DEL PAPER ====\n")

#1a. Ver si el margen de error reportado (RMSEA) coincide con los estadísticos de Chi-cuadrado publicados.
# N/2 = 1178 y con los grados de libertad (cargas + errores + correlaciones entre factores)

#calcula grados de libertad, num_preguntas = p y num_factores = k
df_cfa <- function(p, k) p * (p + 1) / 2 - (2 * p + ifelse(k > 1, k * (k - 1) / 2, 0))
#tabla con los valores reportados en el paper
aud <- data.frame(modelo = c("M1 3F-18", "M2 1F-18", "M3 3F-17"),
                  p = c(18, 18, 17), k = c(3, 1, 3),
                  chisq = c(536.41, 2607.54, 495.91),
                  rmsea_pub = c(.051, .125, .053))
#aplica la función para ver si hay coherencia con los grados de libertad que dijeron
aud$df         <- mapply(df_cfa, aud$p, aud$k)
aud$rmsea_calc <- round(sqrt(pmax(aud$chisq - aud$df, 0) / (aud$df * (N / 2 - 1))), 3)
print(aud)


# 1b. "X" de Mann-Whitney es Z; r de Rosenthal = |Z| / sqrt(N)
z <- c(sexo_tot = 4.344, sexo_rech = 1.762, sexo_apoy = 9.236, sexo_pred = 4.884,
       niv_tot = 13.704, niv_rech = 11.707, niv_apoy = 12.228, niv_pred = 8.519)
# Convierte esos valores a $r$ (un indicador universal de qué tan "fuerte" 
# es la diferencia, independientemente de la muestra.
cat("\nr = Z/sqrt(N):\n"); print(round(z / sqrt(N), 2))
cat("(publicado: .08 .03 .19 .10 | .28 .24 .25 .17)\n")



# 1c. Puntaje total implícito por las medias de factor (Anexo 1: 9 ítems
#     invertidos de Rechazo + 5 Apoyo + 3 Predisposición) vs. publicado
#invierte la media de los items de rechazo con 6 - respuesta, se repite para 5 grupos
tot_imp <- sapply(1:5, function(j) (9 * (6 - tg[[1]]$mu[j]) + 5 * tg[[2]]$mu[j] + 3 * tg[[3]]$mu[j]) / 17)
chk <- data.frame(grupo = names(tot_pub), total_pub = tot_pub, total_implicado = round(tot_imp, 3))
chk$dif <- round(chk$total_pub - chk$total_implicado, 3)
chk$media_item19_si_total_de_18 <- round(18 * chk$total_pub - 17 * tot_imp, 2)
cat("\nTotal publicado vs implicado por las medias de factor:\n"); print(chk)
cat("Dif. sistemática ~ +.04: HIPÓTESIS (no verificable): el total publicado pudo\n",
    "calcularse con 18 ítems (incluye el 19)'.\n")

# 1d. Regla de depuración del paper (|carga|>=.40 y brecha >=.10 vs 2a carga)
# la 2da regla la pusieron para evitar que 1 item explique 2 factores
prune_rule <- function(L, cut = .40, gap = .10) {
  o <- t(apply(abs(L), 1, sort, decreasing = TRUE))
  o[, 1] >= cut & (o[, 1] - o[, 2]) >= gap
}
kept <- tab1$item[prune_rule(Lam)]
cat("\nLa regla reproduce exactamente los 18 ítems del AFE: ",
    setequal(kept, ITEMS18), " | eliminados: ", paste(sort(setdiff(1:24, kept)), collapse = ","), "\n")




## ---------------------------------------------------------------------------
## 2. MODELO POBLACIONAL IMPLÍCITO EN LA TABLA 1  (A1)
## Evalua comunalidades
##        confiabilidad (Alfa de Cronbach) 
##        El ajuste de los modelos confirmatorios (AFC).
## ---------------------------------------------------------------------------
cat("\n==== 2. MODELO POBLACIONAL IMPLÍCITO ====\n")
#En el modelo factorial clásico con factores ortogonales la correlación entre dos preguntas
# $i$ y $j$ se explica únicamente por la suma del producto de sus cargas factoriales

#$$R = \Lambda \Lambda^\transpose + \Theta$$
#Calcula las correlaciones de los 24 items
R_pop <- Lam %*% t(Lam); diag(R_pop) <- 1
dimnames(R_pop) <- list(nm(1:24), nm(1:24))
#Comunalidad a partir de las cargas factoriales
h2 <- rowSums(Lam^2)
#Solo es posible si la matriz R_pop es definida positiva (todo eigenvalor > 0)
#Si hubiera negativo diría que hay datos menos dispersos que en 1 punto, imposible
cat("Comunalidad máx:", round(max(h2), 3), "| menor autovalor de R:", round(min(eigen(R_pop)$values), 3), "\n")

alpha_std <- function(S) { k <- nrow(S); rb <- (sum(S) - k) / (k * (k - 1)); k * rb / (1 + (k - 1) * rb) }
sgn <- c(rep(-1, 9), rep(1, 8))                       # Rechazo se invierte en el total
S17 <- R_pop[nm(ITEMS17), nm(ITEMS17)] * outer(sgn, sgn)
al <- c(F1 = alpha_std(R_pop[nm(F1), nm(F1)]), F2 = alpha_std(R_pop[nm(F2), nm(F2)]),
        F3 = alpha_std(R_pop[nm(F3), nm(F3)]), Total = alpha_std(S17))
cat("Alfas implícitos por la Tabla 1:", round(al, 3), " | publicados: .861 .736 .704 .856\n")
cat("Varianza de los 4 primeros autovalores (%):",
    round(100 * sum(eigen(R_pop)$values[1:4]) / 24, 1), " | publicado: 47.97\n")

mod_txt <- function(f1, f2, f3)
  paste0("RECH =~ ", paste(nm(f1), collapse = " + "), "\n",
         "APOY =~ ", paste(nm(f2), collapse = " + "), "\n",
         "PRED =~ ", paste(nm(f3), collapse = " + "))
MOD <- list(M1 = mod_txt(F1, c(F2, 19), F3),
            M2 = paste0("ATN =~ ", paste(nm(ITEMS18), collapse = " + ")),
            M3 = mod_txt(F1, F2, F3))
vars_of <- function(m) if (m == "M3") ITEMS17 else ITEMS18
for (m in names(MOD)) {
  S <- R_pop[nm(vars_of(m)), nm(vars_of(m))]
  f <- cfa(MOD[[m]], sample.cov = S, sample.nobs = N / 2, std.lv = TRUE)
  fm <- fm <- fitMeasures(f, c("chisq", "df", "gfi", "agfi", "rmr", "rmsea")); Fml <- unname(fm["chisq"]) / (N / 2 - 1)
  cat(sprintf("%s  RMSEA poblacional (ML) = %.3f", m, sqrt(Fml / unname(fm["df"]))))
  if (m != "M2") {
    ss <- standardizedSolution(f)
    ph <- ss[ss$op == "~~" & ss$lhs != ss$rhs & ss$lhs %in% c("RECH", "APOY", "PRED"), ]
    cat("  | corr. entre factores implícitas:",
        paste0(ph$lhs, "-", ph$rhs, "=", round(ph$est.std, 2), collapse = "; "))
  }
  cat("\n")
}
cat("(Publicado, ADF muestral: .051 / .125 / .053.)\n")

## ---------------------------------------------------------------------------
## 3. CALIBRACIÓN DE MARGINALES Y DIFERENCIAS POR GRUPO (forma cerrada)
##    ítem = clip(round(m + s * z)), z ~ N(0,1);  m = m0 + dS*hombre + dE*bach
## ---------------------------------------------------------------------------
cat("\n==== 3. CALIBRACIÓN ====\n")
#La actitud hacia algo es continua, al obligar a tomar un item tipo Lickert
#el alumno discretizó su respuesta, por lo que aqui se corrige
#si tuviese una opinion entre 2.5 y 3.5, marca 3, etc
cuts   <- c(-Inf, 1.5, 2.5, 3.5, 4.5, Inf)
#probabilidad de que una persona elija cada una de las 5 opciones 
#si su actitud latente sigue una distribución Normal con media $m$ y desviación $s$:
cell_p <- function(m, s) diff(pnorm((cuts - m) / s))
#Resta las probabilidades acumuladas consecutivas para obtener el área bajo la curva
#de cada categoría. Es decir: $[P(X=1), P(X=2), P(X=3), P(X=4), P(X=5)]$.
mom_p  <- function(p) { mu <- sum(1:5 * p); c(mu, sqrt(sum((1:5)^2 * p) - mu^2)) }

#Como asumimos independencia en genero y nivel educativo, crea las 4 opciones
#H_bach,H_lic,M_bach,M_licS
cells <- expand.grid(h = c(0, 1), b = c(0, 1))
cells$w <- ifelse(cells$h == 1, p_hombre, 1 - p_hombre) * ifelse(cells$b == 1, p_bach, 1 - p_bach)

#Simula las medias y SD de los 5 grupos a partir de 4 parámetros latentes:
grp_moments <- function(th) { 
  #Calcula la probabilidad de responder de 1 a 5 para cada una de las 4 celdas demográficas:
  P <- t(mapply(function(h, b) cell_p(th[1] + th[3] * h + th[4] * b, exp(th[2])),
                cells$h, cells$b))
  
  #Pondera las probabilidades del subgrupo seleccionado y extrae su media y SD observables:
  agg <- function(sel) { 
    w <- cells$w * sel
    mom_p(colSums(P * w) / sum(w)) 
  }
  
  #Apila los resultados simulados para General, Hombres, Mujeres, Bachillerato y Licenciatura:
  rbind(gen = agg(rep(1, 4)), 
        H   = agg(cells$h == 1), 
        M   = agg(cells$h == 0),
        B   = agg(cells$b == 1), 
        L   = agg(cells$b == 0))
}

#Desviación teórica del promedio del factor si cada ítem tuviera varianza 1:
#Permite corregir el hecho de que promediar varios ítems siempre reduce la dispersión observada.
unit_sd <- sapply(list(F1, F2, F3), function(f) sqrt(sum(R_pop[nm(f), nm(f)])) / length(f))

#Función que busca y calibra los parámetros latentes óptimos para un factor f:
fit_factor <- function(f) {
  mu_t <- tg[[f]]$mu                             #Vector de medias reportadas en el paper para los 5 grupos
  sd_t <- tg[[f]]$sd / unit_sd[f]                #SD objetivo escalada para un ítem individual
  
  #Función de pérdida: penaliza al cuadrado la distancia entre lo simulado y lo publicado:
  obj <- function(th) { 
    g <- grp_moments(th)
    sum(((g[, 1] - mu_t) / .01)^2) + sum(((g[, 2] - sd_t) / .02)^2) 
  }
  
  best <- NULL
  #Búsqueda desde 8 puntos iniciales distintos para evitar quedar atrapado en mínimos locales:
  for (s0 in c(1, 1.5, 2, 3)) for (dm in c(0, -2)) {
    #Ajusta th = (m0, log_s, dS, dE) hasta minimizar el error en obj:
    r <- optim(c(mu_t[1] + dm, log(s0), 0, 0), obj, method = "Nelder-Mead",
               control = list(maxit = 5000, reltol = 1e-10))
    #Conserva la solución con la menor discrepancia:
    if (is.null(best) || r$value < best$value) best <- r
  }
  best
}

#Ejecuta el proceso de calibración para los 3 factores (Rechazo, Apoyo y Predisposición):
fits  <- lapply(1:3, fit_factor)

#Extrae los parámetros calibrados y los organiza en una tabla con nombres descriptivos:
theta <- t(sapply(fits, function(r) r$par))
colnames(theta) <- c("m0", "log_s", "dS", "dE")
rownames(theta) <- c("Rechazo", "Apoyo", "Predisp")

#Muestra los parámetros obtenidos y el error final (valores menores a 20 indican excelente ajuste):
cat("Parámetros calibrados:\n"); print(round(theta, 3))
cat("Valor de la función objetivo (esperable < ~20):", round(sapply(fits, function(r) r$value), 1), "\n")

# Asigna a cada uno de los 24 ítems los parámetros correspondientes según la familia a la que pertenece:
P <- data.frame(m0 = theta[fam, 1], s = exp(theta[fam, 2]), dS = theta[fam, 3], dE = theta[fam, 4])

## ---------------------------------------------------------------------------
## 4. CÓPULA: elegir la correlación latente para que la Pearson ORDINAL = R_pop
##    (números aleatorios comunes => convergencia estable)
## ---------------------------------------------------------------------------
##Calibra iterativamente la correlación latente mediante una cópula gaussiana
##para que, al redondear y recortar a la escala Likert (1-5), 
##las correlaciones ordinales observadas reproduzcan con exactitud la matriz de Pearson
##publicada (R_pop), compensando la pérdida de correlación por discretización.
sim_ord <- function(Rlat, Zb, h, b, P) {
  n   <- nrow(Zb)
  Zl  <- Zb %*% chol(Rlat)
  loc <- matrix(P$m0, n, 24, byrow = TRUE) + outer(h, P$dS) + outer(b, P$dE)
  X   <- round(loc + Zl * matrix(P$s, n, 24, byrow = TRUE))
  X[X < 1] <- 1; X[X > 5] <- 5
  colnames(X) <- nm(1:24); X
}
near_corr <- function(A) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  B <- e$vectors %*% diag(pmax(e$values, 1e-3)) %*% t(e$vectors)
  d <- sqrt(diag(B)); B / outer(d, d)
}
nC <- 100000
Zb <- matrix(rnorm(nC * 24), nC, 24)
hC <- rbinom(nC, 1, p_hombre); bC <- rbinom(nC, 1, p_bach)
Rlat <- R_pop
for (it in 1:12) {
  Ro  <- cor(sim_ord(Rlat, Zb, hC, bC, P))
  err <- max(abs(Ro - R_pop)[upper.tri(Ro)])
  cat(sprintf("  iter %2d  max|R_ordinal - R_Tabla1| = %.4f\n", it, err))
  Rlat <- near_corr(Rlat + 0.8 * (R_pop - Ro))
}
Ro <- cor(sim_ord(Rlat, Zb, hC, bC, P))
cat(sprintf("  final    max|R_ordinal - R_Tabla1| = %.4f (esperable < .002)\n",
            max(abs(Ro - R_pop)[upper.tri(Ro)])))




## ---------------------------------------------------------------------------
## 5. GENERAR UNA MUESTRA SINTÉTICA DE N = 2356
## valida que la muestra reproduzca fielmente las medias, desviaciones, 
## categorías diagnósticas, alfas de Cronbach y tamaños del efecto publicados en el artículo,
## exportando la base a un archivo CSV.
## ---------------------------------------------------------------------------
gen_data <- function() {
  n_h <- round(p_hombre * N)
  h <- sample(c(rep(1, n_h), rep(0, N - n_h)))
  b <- sample(c(rep(1, n_bach), rep(0, n_lic)))
  X <- sim_ord(Rlat, matrix(rnorm(N * 24), N, 24), h, b, P)
  d <- as.data.frame(X)
  d$sexo  <- factor(ifelse(h == 1, "H", "M"))
  d$nivel <- factor(ifelse(b == 1, "Bach", "Lic"))
  d
}
score <- function(d) {
  data.frame(rech = rowMeans(d[, nm(F1)]), apoy = rowMeans(d[, nm(F2)]), pred = rowMeans(d[, nm(F3)]),
             tot  = rowMeans(cbind(6 - d[, nm(F1)], d[, nm(F2)], d[, nm(F3)])))
}
alpha_raw <- function(X) { X <- as.matrix(X); k <- ncol(X); C <- cov(X); k / (k - 1) * (1 - sum(diag(C)) / sum(C)) }
cat3 <- function(x) round(100 * c(mean(x < 2.6), mean(x >= 2.6 & x < 3.4), mean(x >= 3.4)), 1)
rosenthal <- function(x, g) {                    # r = |Z| / sqrt(N) con Z obtenido del p de Mann-Whitney
  p <- suppressWarnings(wilcox.test(x ~ g)$p.value)
  abs(qnorm(p / 2)) / sqrt(length(x))
}

dat <- gen_data()
sc  <- score(dat)

cat("\n==== 5. MUESTRA SINTÉTICA vs. TABLAS PUBLICADAS ====\n")
t3 <- data.frame(medida = c("Total", "Rechazo", "Apoyo", "Predisp"),
                 M_sim  = round(colMeans(sc[, c("tot", "rech", "apoy", "pred")]), 2),
                 M_pub  = c(1.98, 3.65, 1.35, 1.70),
                 SD_sim = round(apply(sc[, c("tot", "rech", "apoy", "pred")], 2, sd), 3),
                 SD_pub = c(.659, .945, .656, .723))
cat("Tabla 3\n"); print(t3, row.names = FALSE)
cat("\nTabla 4 (% bajo/medio/alto):\n")
print(rbind(Rechazo_sim = cat3(sc$rech), Rechazo_pub = c(13.6, 23.3, 63.1),
            Apoyo_sim = cat3(sc$apoy),   Apoyo_pub = c(92.4, 5.7, 1.9),
            Predisp_sim = cat3(sc$pred), Predisp_pub = c(84.8, 11.6, 3.5)))
cat("\nAlfas (sim): F1 =", round(alpha_raw(dat[, nm(F1)]), 3), " F2 =", round(alpha_raw(dat[, nm(F2)]), 3),
    " F3 =", round(alpha_raw(dat[, nm(F3)]), 3),
    " Total =", round(alpha_raw(cbind(6 - dat[, nm(F1)], dat[, nm(F2)], dat[, nm(F3)])), 3),
    "\n       (pub): .861 .736 .704 .856\n")
cat("r Rosenthal sexo  (tot,rech,apoy,pred):", round(sapply(sc[, c("tot","rech","apoy","pred")], rosenthal, g = dat$sexo), 2), " pub: .08 .03 .19 .10\n")
cat("r Rosenthal nivel (tot,rech,apoy,pred):", round(sapply(sc[, c("tot","rech","apoy","pred")], rosenthal, g = dat$nivel), 2), " pub: .28 .24 .25 .17\n")

write.csv(dat, "ATN_sintetico.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 6. RE-EJECUCIÓN DEL PIPELINE DEL PAPER
##    mitad 1 -> AFE (PAF, Equamax, Pearson) ; mitad 2 -> AFC (3 modelos)
## ---------------------------------------------------------------------------
align_to <- function(L, Lref) {                  # permuta/invierte signo para comparar con Tabla 1
  cg <- psych::factor.congruence(L, Lref)
  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  best <- perms[[which.max(sapply(perms, function(p) sum(abs(diag(cg[p, ])))))]]
  sweep(L[, best, drop = FALSE], 2, sign(diag(cg[best, ])), "*")
}
cfa_fit <- function(model, data, est) {
  f <- suppressWarnings(tryCatch(cfa(model, data = data, estimator = est, std.lv = TRUE),
                                 error = function(e) NULL))
  out <- setNames(rep(NA_real_, 5), c("chisq", "gfi", "agfi", "rmr", "rmsea"))
  if (is.null(f) || !lavInspect(f, "converged")) return(out)
  fm <- fitMeasures(f, c("chisq", "gfi", "agfi", "rmr", "rmsea"))
  out[] <- unname(fm[c("chisq", "gfi", "agfi", "rmr", "rmsea")])   # NA si la medida no existe
  out
}

run_pipeline <- function(d, do_pa = TRUE, ests = c("WLS", "ML")) {
  ord <- sample(nrow(d)); n1 <- nrow(d) %/% 2
  h1 <- d[ord[1:n1], ]; h2 <- d[ord[(n1 + 1):nrow(d)], ]
  X1 <- h1[, nm(1:24)]; R1 <- cor(X1)

  ev    <- eigen(R1, symmetric = TRUE)$values
  efa   <- psych::fa(X1, nfactors = 3, fm = "pa", rotate = "equamax")
  La    <- align_to(unclass(efa$loadings)[nm(1:24), ], Lam)
  kept  <- which(prune_rule(La))
  nfpa  <- NA
  if (do_pa) { invisible(capture.output(pa <- psych::fa.parallel(X1, fa = "fa", n.iter = 30, plot = FALSE)))
               nfpa <- pa$nfact }

  s <- c(KMO = psych::KMO(R1)$MSA, bartlett_p = psych::cortest.bartlett(R1, n = n1)$p.value,
         n_kaiser = sum(ev > 1), n_PA = nfpa, var4 = 100 * sum(ev[1:4]) / 24,
         rmse_cargas = sqrt(mean((La - Lam)^2)), mismos_18 = as.numeric(setequal(kept, ITEMS18)),
         alpha_F1 = alpha_raw(d[, nm(F1)]), alpha_F2 = alpha_raw(d[, nm(F2)]), alpha_F3 = alpha_raw(d[, nm(F3)]),
         alpha_tot = alpha_raw(cbind(6 - d[, nm(F1)], d[, nm(F2)], d[, nm(F3)])))
  for (e in ests) for (m in names(MOD)) {
    r <- cfa_fit(MOD[[m]], h2[, nm(vars_of(m))], e)
    names(r) <- paste(e, m, names(r), sep = "_"); s <- c(s, r)
  }
  list(summary = s, efa = efa, loadings = La, half2 = h2)
}

cat("\n==== 6. PIPELINE SOBRE LA MUESTRA SINTÉTICA ====\n")
res <- run_pipeline(dat)
print(round(res$summary, 3))
cat("\nCargas recuperadas (alineadas) vs Tabla 1, ítems de la escala final:\n")
print(round(cbind(sim = res$loadings[nm(ITEMS18), ], pub = Lam[nm(ITEMS18), ]), 2))

# Alternativa moderna: ítems tratados como ordinales (correlaciones policóricas, WLSMV)
f_ord <- tryCatch(cfa(MOD$M3, data = res$half2[, nm(ITEMS17)], ordered = nm(ITEMS17),
                      estimator = "WLSMV", std.lv = TRUE), error = function(e) NULL)
if (!is.null(f_ord)) {
  fm <- fitMeasures(f_ord)
  cat("\nM3 con WLSMV (ordinal):\n")
  print(round(fm[intersect(c("chisq.scaled", "df.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"),
                           names(fm))], 3))
  ss <- standardizedSolution(f_ord)
  cat("Correlaciones entre factores (WLSMV):\n")
  print(ss[ss$op == "~~" & ss$lhs != ss$rhs & ss$lhs %in% c("RECH", "APOY", "PRED"),
           c("lhs", "rhs", "est.std", "ci.lower", "ci.upper")])
}

## ---------------------------------------------------------------------------
## 7. MONTE CARLO: ¿están los valores publicados dentro de lo plausible?
##    (chequeo predictivo: B muestras sintéticas completas, nuevo split cada vez)
## ---------------------------------------------------------------------------
B <- 100                                          # súbelo a 500+ si hay tiempo
cat("\n==== 7. MONTE CARLO (B =", B, ") ====\n")
mc <- t(replicate(B, run_pipeline(gen_data(), do_pa = FALSE)$summary))

pub <- c(KMO = .906, var4 = 47.97, alpha_F1 = .861, alpha_F2 = .736, alpha_F3 = .704, alpha_tot = .856,
         WLS_M1_chisq = 536.41, WLS_M1_gfi = .932, WLS_M1_agfi = .912, WLS_M1_rmr = .129, WLS_M1_rmsea = .051,
         WLS_M2_chisq = 2607.54, WLS_M2_gfi = .707, WLS_M2_agfi = .629, WLS_M2_rmr = .145, WLS_M2_rmsea = .125,
         WLS_M3_chisq = 495.91, WLS_M3_gfi = .934, WLS_M3_agfi = .913, WLS_M3_rmr = .122, WLS_M3_rmsea = .053)
tab_mc <- t(sapply(names(pub), function(k) {
  x <- mc[, k]; q <- quantile(x, c(.025, .5, .975), na.rm = TRUE)
  c(publicado = unname(pub[k]), mc_p2.5 = q[[1]], mc_mediana = q[[2]], mc_p97.5 = q[[3]],
    dentro = as.numeric(pub[k] >= q[[1]] & pub[k] <= q[[3]]), n_NA = sum(is.na(x)))
}))
print(round(tab_mc, 3))
cat("\nAFE: % de réplicas con regla de depuración = 18 ítems del paper:", round(100 * mean(mc[, "mismos_18"]), 1), "\n")
cat("AFE: % de réplicas con 4 autovalores > 1:", round(100 * mean(mc[, "n_kaiser"] >= 4), 1),
    "| RMSE medio de cargas vs Tabla 1:", round(mean(mc[, "rmse_cargas"]), 3), "\n")