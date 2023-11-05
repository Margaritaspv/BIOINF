
# Instalar RCurl
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Antes de cargar los datos, comprobar que estamos en el directorio correcto!
getwd()

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/datos-trabajoR", followlocation = TRUE)
data = read.table(file = "datos-trabajoR.txt", head = T)
data #Para comprobar que la tabla es correcta 

# Emplea las funciones head(), summary(), dim() y str(). 
#¿Cuántas variables hay? ¿Cuántos tratamientos? --> Hay 2 variables, y 5 tratamientos.
head (data)
summary (data)
dim(data)
str(data)

#Haz un boxplot para nuestros datos. Uno para cada variable. Colorea Variable 1 y a Variable 2 de forma diferente.
boxplot(Variable1~Tratamiento, data=data, main = "BloxPlot Variable 1", col = "violet")
boxplot(Variable2~Tratamiento, data=data, main = "BloxPlot Variable 2", col = "yellow")

#Haz un gráfico de dispersion con las dos variables. Cada tratamiento debe ir de un color distinto. 
plot (x = data$Variable1, y = data$Variable2, main = "Scatter plot de Variables", col = data$Tratamiento, xlab = "Variable 1", ylab = "Variable 2")

#Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho.
legend ("bottomright", legend = c("Tratamiento 1", "Tratamiento 2", "Tratamiento 3", "Tratamiento 4", "Tratamiento 5"), col = c("black","red","green","turquoise","blue"), pch = 16)

#Haz un histograma para cada variable.
hist (data$Variable1, col = "violet", main = "Histograma Variable 1", xlab = "Variable 1")
hist (data$Variable2, col = "yellow", main = "Histograma Variable 2", xlab= "Variable 2")

#Haz un factor en la columna de tratamiento y guárdalo en una variable
factor_tratamiento <- factor(data$Tratamiento)
factor_tratamiento

#Calcula la media y la desviación estándar para cada tratamiento.
mean1 <- aggregate (Variable1 ~ factor_tratamiento, data, mean)
mean1
mean2 <- aggregate (Variable2 ~factor_tratamiento, data, mean)
mean2
	
       #Esto me daría de resultado la media y la desviación de cada variable
media_desviacion <- aggregate(. ~ Tratamiento, data, 
                               FUN = function(x) c(media = mean(x), desviacion = sd(x)))
media_desviacion

	 #Esto me daría de resultado la media y la desviación de ambas variables
mean_variable1 = tapply(data$Variable1, data$Tratamiento, mean)
mean_variable2 = tapply(data$Variable2, data$Tratamiento, mean)
mean_total = (mean_variable1 + mean_variable2) / 2
mean_total
sd_variable1 = tapply(data$Variable1, data$Tratamiento, sd)
sd_variable2 = tapply(data$Variable2, data$Tratamiento, sd)
sd_total = sqrt(sd_variable1^2 + sd_variable2^2)
sd_total

#Averigua cuántos elementos tiene cada tratamiento.
elementos = table(data$Tratamiento)
elementos

#Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
bp = read.table("datos-trabajoR.txt",header=T)
Tratamiento1 = data[data$Tratamiento == 1,]
Tratamiento1
Tratamiento4 = data[data$Tratamiento == 4,]
Tratamiento4

#Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la
#Variable 1 son iguales. ¿Puedes comprobarlo? Para ello, necesitarás comprobar
#primero si los datos se distribuyen de forma normal. En función del resultado de la
#prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras
#son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus
#resultados.

shapiro.test(Tratamiento1$Variable1)
shapiro.test(Tratamiento4$Variable1)
#Los pvalue no son estadisticamente significativos, y siguen distribuciónr normal. Por eso usamos t.test.


t.test(Tratamiento1$Variable1, Tratamiento4$Variable1)
#Los resultados de las medias son diferentes.

var.test(Tratamiento1$Variable1, Tratamiento4$Variable1) 
#El pvalue pequeño y el intervalo de confianza nos indican que las varianzas son diferentes. El ratio tan pequeño y diferente de 1, nos lo corrobora.