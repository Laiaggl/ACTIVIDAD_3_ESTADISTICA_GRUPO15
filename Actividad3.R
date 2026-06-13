rm(list = ls())

set.seed(123)

styleColor = "#00cec9"
p_significance_threshold = 0.05;

paquetes <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "factoextra",
  "pheatmap",
  "gtsummary",
  "broom",
  "pROC"
)

for (p in paquetes) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}
### ============================================================================
### 2. CARGA, ESTRUCTURA Y LIMPIEZA DE LOS DATOS
### ============================================================================


datos = read.csv("Base de datos Los Simpson completa.csv")

genes = c("ADCY3","AGRP","ANKRD27","ANO4","BDNF","CADM1","CADM2","CALCR","CNR1","CREBRF","DPP9","FTO","GIPR","GPR151","GPR75","KIAA0586","KIAA1109","KSR2","LEP","LEPR","MC4R","MRAP2","NRP1","NRP2","NTRK2","PCSK1","PDE3B","PHIP","POMC","PPARG","ROBO1","SH2B1","SIM1","SPARC","TMEM18","UBR2","UHMK1")

## just making sure that the columns that the char columns are trimmed and normalized with Capital Case
datos <- datos %>%
  mutate(across(where(is.character), ~ str_to_title(str_trim(.x))))

datos <- datos %>%
  mutate(
    sexo = factor(
      sexo,
      levels = unique(datos$sexo)
    ),
    
    tabaco = factor(
      tabaco,
      levels = unique(datos$tabaco)
    ),
    
    alcohol = factor(
      alcohol,
      levels = unique(datos$alcohol)
    ),
    
    estado_civil = factor(
      estado_civil,
      levels = unique(datos$estado_civil)
    ),
    
    enfermedades = factor(
      enfermedades,
      levels = unique(datos$enfermedades)
    ),
    
    medicamento = factor(
      medicamento,
      levels = unique(datos$medicamento)
    ),
    
    ansiedad = factor(
      ansiedad,
      levels = unique(datos$ansiedad)
    ),
    
    depresion = factor(
      depresion,
      levels = unique(datos$depresion)
    ),
    
    trastornos = factor(
      trastornos,
      levels = unique(datos$trastornos)
    ),
    
    sobrepeso = factor(
      case_when(
      imc_kg_m2 < 25  ~ "No",
      imc_kg_m2 >= 25 ~ "Sí"
    ),
    levels = c("No", "Sí")
    )
  )

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### Antes de aplicar un PCA se debe asegurar que:
### - Solo se incluyen las variables de expresion seleccionadas.
### - No quedan valores perdidos en la matriz utilizada.
### - Todas las variables se encuentran en formato numerico.
### - El numero de registros conservados tras la limpieza es razonable.
### ============================================================================

colSums(
  is.na(
    datos[, genes]
  )
)

n_inicial <- nrow(datos)

datos_limpios <- datos %>%
  drop_na(
    all_of(genes)
  )

n_final <- nrow(datos_limpios)

data.frame(
  registros_iniciales = n_inicial,
  registros_finales = n_final,
  registros_eliminados = n_inicial - n_final
)

any(
  is.na(
    datos_limpios[, genes]
  )
)

summary(
  datos_limpios[, c(
    "edad_anios",
    "imc_kg_m2",
    "ccintura_cm"
  )]
)


### ============================================================================
### 3. COMPROBACION DE NORMALIDAD EN LA EXPRESION GENICA
### ============================================================================

### El test de Shapiro-Wilk evalua si los datos muestran evidencia contra
### una distribucion normal.
###
### Hipotesis nula: la variable sigue una distribucion normal.
### Si p < 0.05, se rechaza la normalidad con el criterio habitual.
### Si p >= 0.05, no se dispone de evidencia suficiente para rechazarla.
###
### No rechazar la normalidad no equivale a demostrar que la distribucion
### sea perfectamente normal.

normalidad_genes <- lapply(
  genes,
  function(gen) {
    
    test <- shapiro.test(
      datos_limpios[[gen]]
    )
    
    data.frame(
      gen = gen,
      W = unname(test$statistic),
      p_valor = test$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    W = round(W, 4),
    p_valor = round(p_valor, 4),
    
    decision = ifelse(
      p_valor < 0.05,
      "Evidencia de desviacion de normalidad",
      "Sin evidencia suficiente contra normalidad"
    )
  )

normalidad_genes

table(normalidad_genes$decision)

write.csv(
  normalidad_genes,
  file = "resultado_normalidad_panel_genico.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### La normalidad es especialmente relevante para decidir como describir
### los valores de expresion:
### - Media y desviacion estandar cuando la distribucion sea razonable.
### - Mediana y rango intercuartil cuando exista asimetria o evidencia clara
###   de no normalidad.
###
### El PCA no exige normalidad perfecta. En cambio, si los genes tienen
### escalas diferentes, centrar y escalar las variables resulta esencial.
### ============================================================================


### ============================================================================
### 4. ANALISIS DE COMPONENTES PRINCIPALES
### ============================================================================

matriz_expresion <- datos_limpios[, genes]

pca <- prcomp(
  matriz_expresion,
  center = TRUE,
  scale. = TRUE
)

summary(pca)

### --------------------------------------------------------------------------
### 4.1. Varianza explicada
### --------------------------------------------------------------------------

varianza <- pca$sdev^2

proporcion <- varianza / sum(varianza)

tabla_varianza <- data.frame(
  componente = paste0(
    "PC",
    seq_along(proporcion)
  ),
  
  varianza_explicada = round(
    proporcion * 100,
    2
  ),
  
  varianza_acumulada = round(
    cumsum(proporcion) * 100,
    2
  )
)


tabla_varianza[1:8, ]

write.csv(
  tabla_varianza,
  file = "pca_varianza_explicada.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

fviz_eig(
  pca,
  addlabels = TRUE
) +
  labs(
    title = "Varianza explicada por los componentes principales",
    x = "Componente principal",
    y = "Varianza explicada (%)"
  )

### --------------------------------------------------------------------------
### 4.2. Scores de los registros
### --------------------------------------------------------------------------

scores <- as.data.frame(
  pca$x[, 1:6, drop = FALSE]
)

datos_pca <- bind_cols(
  datos_limpios,
  scores
)

head(
  datos_pca[, c(
    "X",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6"
  )]
)

### --------------------------------------------------------------------------
### 4.3. Cargas de los genes
### --------------------------------------------------------------------------

cargas <- as.data.frame(
  pca$rotation[, 1:6, drop = FALSE]
) %>%
  tibble::rownames_to_column(
    var = "gen"
  )

cargas

cargas_PC1 <- cargas %>%
  mutate(
    carga_absoluta = abs(PC1)
  ) %>%
  arrange(
    desc(carga_absoluta)
  )

cargas_PC2 <- cargas %>%
  mutate(
    carga_absoluta = abs(PC2)
  ) %>%
  arrange(
    desc(carga_absoluta)
  )

cargas_PC1

cargas_PC2

write.csv(
  cargas,
  file = "pca_cargas_genes.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### Las cargas permiten interpretar el significado del componente.
###
### Una carga positiva elevada indica que una mayor expresion del gen
### incrementa el score del componente.
###
### Una carga negativa elevada en valor absoluto indica que una mayor
### expresion del gen se relaciona con scores mas bajos.
###
### La interpretacion biologica de PC1 o PC2 debe apoyarse en los genes
### con mayor carga absoluta, no solo en el porcentaje de varianza explicado.
### ============================================================================


### ============================================================================
### 5. REPRESENTACION DEL PCA
### ============================================================================

### --------------------------------------------------------------------------
### 5.1. Variables en el plano PC1-PC2
### --------------------------------------------------------------------------

fviz_pca_var(
  pca,
  col.var = "cos2",
  gradient.cols = c("#2E86AB", "#F6C85F", "#D1495B"),
  repel = TRUE
) +
  labs(
    title = "Representacion de genes en el plano PC1-PC2",
    color = "Cos2"
  )

### El cos2 indica la calidad de representacion del gen en el plano mostrado.
### Genes con cos2 elevado quedan bien representados mediante PC1 y PC2.
### Genes con cos2 bajo pueden estar mejor explicados por otras componentes.


### --------------------------------------------------------------------------
### 5.2. Contribucion de los genes
### --------------------------------------------------------------------------

fviz_contrib(
  pca,
  choice = "var",
  axes = 1,
  top = 10
) +
  labs(
    title = "Genes con mayor contribucion a PC1"
  )

fviz_contrib(
  pca,
  choice = "var",
  axes = 2,
  top = 10
) +
  labs(
    title = "Genes con mayor contribucion a PC2"
  )


### --------------------------------------------------------------------------
### 5.3. Individuos coloreados segun alteracion glucemica
### --------------------------------------------------------------------------

datos_pca <- datos_pca %>%
  mutate(
    categoria_imc = case_when(
      imc_kg_m2 < 25           ~ "Normal",
      imc_kg_m2 >= 25 & imc_kg_m2 < 30 ~ "Sobrepeso",
      imc_kg_m2 >= 30          ~ "Obesidad"
    )
  )

fviz_pca_ind(
  pca,
  geom.ind = "point",
  habillage = datos_pca$sobrepeso,
  addEllipses = TRUE,
  ellipse.type = "confidence",
  palette = c("#2E86AB", "#F6C85F", "#D1495B"),
  legend.title = "Alteracion glucemica"
) +
  labs(
    title = "Distribucion de registros segun perfil transcriptomico"
  )

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### En el grafico de individuos interesa observar:
### - Si los grupos se separan o se solapan.
### - Si la separacion ocurre principalmente en PC1 o en PC2.
### - Si existen registros alejados del conjunto principal.
###
### La separacion visual puede sugerir diferencias, pero no demuestra
### por si sola una asociacion estadistica.
### ============================================================================


### ============================================================================
### 6. CLUSTERING JERARQUICO
### ============================================================================

### En este ejemplo se utiliza clustering jerarquico y cuatro grupos.
### La distancia se calcula a partir de los scores de las primeras cuatro
### componentes, que resumen la informacion principal de expresion.

scores_cluster <- datos_pca %>%
  select(
    PC1,
    PC2,
    PC3,
    PC4
  ) %>%
  scale()

distancias <- dist(
  scores_cluster,
  method = "euclidean"
)

cluster_hierarquico <- hclust(
  distancias,
  method = "ward.D2"
)

plot(
  cluster_hierarquico,
  labels = FALSE,
  hang = -1,
  main = "Dendrograma de perfiles transcriptomicos",
  xlab = "Registros",
  ylab = "Distancia"
)

rect.hclust(
  cluster_hierarquico,
  k = 4,
  border = 2:5
)

datos_pca$cluster <- factor(
  cutree(
    cluster_hierarquico,
    k = 4
  )
)

table(datos_pca$cluster)

table(
  datos_pca$cluster,
  datos_pca$sobrepeso
)

resumen_cluster <- datos_pca %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    edad_mediana = median(edad_anios),
    imc_mediano = median(imc_kg_m2),
    alteracion_n = sum(sobrepeso == "Si"),
    alteracion_pct = mean(sobrepeso == "Si") * 100,
    .groups = "drop"
  )

resumen_cluster

ggplot(
  datos_pca,
  aes(
    x = PC1,
    y = PC2,
    color = cluster,
    shape = sobrepeso
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "Clusters transcriptomicos en el plano PC1-PC2",
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    shape = "Sobrepeso"
  ) +
  theme_classic()

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### Un cluster identifica similitud en el patron conjunto de expresion.
### No representa automaticamente un diagnostico ni un mecanismo causal.
###
### Tras identificar grupos, resulta necesario describirlos mediante
### variables clinicas para conocer si muestran perfiles diferenciados.
### ============================================================================


### ============================================================================
### 7. HEATMAP DE CORRELACIONES Y EXPRESION
### ============================================================================

### Se explora la relacion entre los genes y las primeras cuatro componentes.
### Se emplea correlacion de Spearman, adecuada para valorar asociaciones
### monotonicas sin asumir normalidad estricta.

correlaciones <- cor(
  datos_pca[, genes],
  datos_pca[, c("PC1", "PC2", "PC3", "PC4")],
  method = "spearman"
)

round(
  correlaciones,
  2
)

pheatmap(
  correlaciones,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  cutree_rows = 3,
  cutree_cols = 3,
  number_format = "%.2f",
  main = "Correlacion entre genes y componentes principales"
)

anotacion <- data.frame(
  sobrepeso = datos_pca$sobrepeso,
  cluster = datos_pca$cluster
)

rownames(anotacion) <- datos_pca$X

colnames(anotacion)


## TODO fix this
## Error in annotation_colors[[colnames(annotation)[i]]] : 
## subscript out of bounds
pheatmap(
  t(
    as.matrix(
      datos_pca[, genes]
    )
  ),
  scale = "row",
  annotation_col = anotacion,
  show_colnames = FALSE,
  main = "Patrones de expresion del panel genico"
)

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### En el heatmap de correlaciones:
### - Valores positivos elevados indican que el gen aumenta con el componente.
### - Valores negativos elevados indican una relacion inversa.
### - Valores proximos a cero indican poca relacion con ese componente.
###
### Los patrones de color permiten reconocer bloques de genes con
### comportamiento semejante, que despues deben interpretarse con cautela.
### ============================================================================


### ============================================================================
### 8. AGRUPACION DE LOS COMPONENTES EN CUARTILES
### ============================================================================

### Para este ejemplo se utilizan cuartiles.
### Un cuartil divide la distribucion en cuatro grupos de tamano aproximado:
### Q1 contiene los valores mas bajos y Q4 los mas altos.
###
### La utilidad de agrupar componentes es facilitar comparaciones descriptivas
### y permitir interpretar modelos mediante categorias de referencia.

crear_cuartiles <- function(x) {
  
  puntos <- quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE
  )
  
  cut(
    x,
    breaks = c(-Inf, puntos, Inf),
    labels = c("Q1", "Q2", "Q3", "Q4"),
    include.lowest = TRUE,
    right = TRUE
  )
}

datos_pca <- datos_pca %>%
  mutate(
    PC1_q = crear_cuartiles(PC1),
    PC2_q = crear_cuartiles(PC2)
  )

table(datos_pca$PC1_q)

table(datos_pca$PC2_q)

datos_pca %>%
  group_by(PC1_q) %>%
  summarise(
    n = n(),
    PC1_min = min(PC1),
    PC1_max = max(PC1),
    PC1_mediana = median(PC1),
    .groups = "drop"
  )

### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### Un valor alto de PC1 no es necesariamente favorable ni desfavorable.
### Su interpretacion depende de las cargas:
### - Si genes inflamatorios tienen cargas positivas elevadas, Q4 puede
###   representar un perfil inflamatorio mas marcado.
### - Si predominan genes protectores con cargas positivas, la lectura puede
###   ser diferente.
### ============================================================================


### ============================================================================
### 9. TABLA DESCRIPTIVA SEGUN CUARTILES DE PC1
### ============================================================================

### Se seleccionan los seis genes con mayor carga absoluta en PC1 para
### construir una tabla compacta e interpretable.

genes_destacados <- cargas_PC1$gen[1:6]

genes_destacados

tabla_genes_pc1 <- datos_pca %>%
  select(
    PC1_q,
    all_of(genes_destacados)
  ) %>%
  tbl_summary(
    by = PC1_q,
    statistic = all_continuous() ~ "{median} ({p25}, {p75})",
    digits = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  add_p(
    test = all_continuous() ~ "kruskal.test",
    pvalue_fun = function(x) {
      style_pvalue(
        x,
        digits = 3
      )
    }
  ) %>%
  modify_caption(
    "**Expresion de genes destacados segun cuartiles de PC1**"
  )

tabla_genes_pc1

tabla_clinica_pc1 <- datos_pca %>%
  select(
    PC1_q,
    edad_anios,
    sexo,
    imc_kg_m2,
    ccintura_cm,
    sobrepeso
  ) %>%
  tbl_summary(
    by = PC1_q,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ 1
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = function(x) {
      style_pvalue(
        x,
        digits = 3
      )
    }
  ) %>%
  modify_caption(
    "**Caracteristicas clinicas segun cuartiles de PC1**"
  )

tabla_clinica_pc1



### ============================================================================
### 10. REGRESION LOGISTICA PARA ALTERACION GLUCEMICA
### ============================================================================

### La regresion logistica se utiliza cuando la variable resultado es binaria.
### En este caso:
### - Evento: sobrepeso = "Si".
### - Referencia: sobrepeso = "No".
###
### Se estudia si los cuartiles de PC1 se asocian con la presencia de
### alteracion glucemica.
###
### Q1 se establece como categoria de referencia.

datos_pca <- datos_pca %>%
  mutate(
    sobrepeso = relevel(
      sobrepeso,
      ref = "No"
    ),
    
    PC1_q = relevel(
      PC1_q,
      ref = "Q1"
    ),
    
    sexo = relevel(
      sexo,
      ref = "Femenino"
    ),
    
    tabaco = relevel(
      tabaco,
      ref = "No Fumador"
    )
  )

modelo_crudo <- glm(
  sobrepeso ~ PC1_q,
  data = datos_pca,
  family = binomial(link = "logit")
)

modelo_demografico <- glm(
  sobrepeso ~ PC1_q + edad_anios + sexo,
  data = datos_pca,
  family = binomial(link = "logit")
)

## TODO review this. it's giving warnings and seems to be incorrect

modelo_ajustado <- glm(
  sobrepeso ~ PC1_q + edad_anios + sexo + tabaco,
  data = datos_pca,
  family = binomial(link = "logit")
)

summary(modelo_crudo)

summary(modelo_demografico)

summary(modelo_ajustado)


### --------------------------------------------------------------------------
### 10.1. Extraccion de OR e intervalos de confianza
### --------------------------------------------------------------------------

extraer_or <- function(modelo, etiqueta) {
  
  broom::tidy(
    modelo,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      modelo = etiqueta,
      
      OR_IC95 = paste0(
        round(estimate, 2),
        " (",
        round(conf.low, 2),
        "; ",
        round(conf.high, 2),
        ")"
      ),
      
      p_valor = ifelse(
        p.value < 0.001,
        "<0.001",
        format(
          round(p.value, 3),
          nsmall = 3
        )
      )
    ) %>%
    select(
      modelo,
      term,
      OR = estimate,
      conf.low,
      conf.high,
      OR_IC95,
      p_valor
    )
}

or_crudo <- extraer_or(
  modelo_crudo,
  "Modelo crudo"
)

or_demografico <- extraer_or(
  modelo_demografico,
  "Ajustado por edad y sexo"
)

or_ajustado <- extraer_or(
  modelo_ajustado,
  "Ajuste ampliado"
)

tabla_or <- bind_rows(
  or_crudo,
  or_demografico,
  or_ajustado
)

tabla_or

tabla_or_pc1 <- tabla_or %>%
  filter(
    grepl(
      "^PC1_q",
      term
    )
  )

tabla_or_pc1

### ============================================================================
### PUNTO DE ATENCION: QUE ES UNA OR
### ============================================================================
### La Odds Ratio compara las odds de presentar el evento entre categorias.
###
### odds = probabilidad del evento / probabilidad de no presentar el evento
###
### OR = 1:
### No se observa diferencia en las odds frente a la categoria de referencia.
###
### OR > 1:
### Se observan mayores odds del evento respecto a la referencia.
###
### OR < 1:
### Se observan menores odds del evento respecto a la referencia.
###
### Ejemplo de redaccion:
### "El cuartil Q4 de PC1 presento mayores odds de alteracion glucemica
### frente a Q1, tras ajustar por edad, sexo, IMC, actividad y tabaquismo."
###
### El intervalo de confianza del 95% debe comprobarse siempre:
### - Si incluye 1, la estimacion es compatible con ausencia de asociacion.
### - Si queda completamente por encima de 1, la asociacion es positiva.
### - Si queda completamente por debajo de 1, la asociacion es inversa.
###
### La OR no debe interpretarse como una diferencia de medias ni como una
### probabilidad directa de presentar el evento.
### ============================================================================


### ============================================================================
### 11. PROBABILIDADES PREDICHAS Y CALIDAD DEL MODELO
### ============================================================================

### Las probabilidades predichas facilitan la interpretacion absoluta del
### modelo, mientras que las OR expresan comparaciones relativas.

perfil_prediccion <- expand.grid(
  PC1_q = factor(
    c("Q1", "Q2", "Q3", "Q4"),
    levels = levels(datos_pca$PC1_q)
  ),
  
  edad_anios = median(datos_pca$edad_anios),
  
  sexo = factor(
    "Femenino",
    levels = levels(datos_pca$sexo)
  ),
  
  imc_kg_m2 = median(datos_pca$imc_kg_m2),
  
  grasa_corporal_pct = median(datos_pca$grasa_corporal_pct),
  
  tabaco = factor(
    "No Fumador",
    levels = levels(datos_pca$tabaco)
  ),
  
  sobrepeso = factor(
    "No",
    levels = levels(datos_pca$sobrepeso)
  )
)

perfil_prediccion$probabilidad <- predict(
  modelo_ajustado,
  newdata = perfil_prediccion,
  type = "response"
)

perfil_prediccion

ggplot(
  perfil_prediccion,
  aes(
    x = PC1_q,
    y = probabilidad
  )
) +
  geom_col() +
  labs(
    title = "Probabilidad estimada de alteracion glucemica",
    subtitle = "Perfil ajustado a valores de referencia",
    x = "Cuartil de PC1",
    y = "Probabilidad estimada"
  ) +
  theme_classic()

roc_modelo <- roc(
  response = datos_pca$sobrepeso,
  predictor = fitted(modelo_ajustado),
  levels = c("No", "Sí"),
  direction = "<",
  quiet = TRUE
)

auc_modelo <- auc(roc_modelo)

auc_modelo

plot(
  roc_modelo,
  main = "Curva ROC del modelo ajustado"
)


### ============================================================================
### PUNTO DE ATENCION
### ============================================================================
### El AUC informa sobre la capacidad discriminativa:
### - AUC de 0.50: capacidad semejante al azar.
### - AUC alrededor de 0.70: discriminacion moderada.
### - AUC alrededor de 0.80 o superior: discriminacion buena.
###
### Un AUC elevado no significa que las probabilidades esten perfectamente
### calibradas. Describe la capacidad para ordenar registros con mayor o
### menor probabilidad del evento.
###
### La conclusion del modelo debe integrar:
### - Sentido biologico del componente.
### - Magnitud de las OR.
### - Intervalos de confianza.
### - Cambio de las estimaciones tras el ajuste.
### - Capacidad discriminativa.
### ============================================================================


### ============================================================================
### 12. TABLA FINAL DEL MODELO AJUSTADO
### ============================================================================

tabla_modelo_final <- tbl_regression(
  modelo_ajustado,
  exponentiate = TRUE,
  label = list(
    PC1_q ~ "Cuartiles de PC1",
    edad_anios ~ "Edad, anos",
    sexo ~ "Sexo",
    tabaco ~ "Tabaquismo"
  )
) %>%
  bold_labels() %>%
  modify_caption(
    "**Modelo ajustado para alteracion glucemica**"
  )

tabla_modelo_final


### ============================================================================
### 13. EXPORTACION DE RESULTADOS
### ============================================================================

write.csv(
  datos_pca,
  file = "base_con_componentes_cuartiles_y_clusters.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  cargas,
  file = "resultado_cargas_componentes.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tabla_or,
  file = "resultado_modelos_logisticos.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  resumen_cluster,
  file = "descripcion_clusters_transcriptomicos.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


### ============================================================================
### 14. ELEMENTOS CENTRALES PARA LA INTERPRETACION
### ============================================================================

### NORMALIDAD
### Revisar que genes presentan evidencia de desviacion de normalidad y
### justificar el tipo de resumen descriptivo seleccionado.

### PCA
### Identificar cuanta varianza explican PC1 y PC2 y que genes tienen las
### cargas absolutas mas elevadas en cada componente.

### REPRESENTACION
### Valorar si los perfiles con alteracion glucemica ocupan regiones
### diferenciadas del plano PC1-PC2 o si existe un elevado solapamiento.

### CLUSTERING
### Describir si los clusters se diferencian en HbA1c, glucosa, IMC o
### frecuencia de alteracion glucemica.

### TABLAS
### Identificar si los valores de expresion y las variables clinicas cambian
### de forma ordenada al pasar de Q1 a Q4 de PC1.

### REGRESION LOGISTICA
### Especificar el evento, la categoria de referencia, las OR, los IC95% y
### las diferencias entre el modelo crudo y el modelo ajustado.

### INTERPRETACION GLOBAL
### No limitar la conclusion al valor p. Integrar la magnitud de los efectos,
### su precision, la plausibilidad biologica y la calidad del modelo.
