library(MaxWiK)
RES = MaxWiK::Results_toy_experiments
RES$method_kernel  =  paste0( RES$method_name, '_', RES$kernel_name )
RES  =  RES[ , 3:9]

RES[ which( RES$method_kernel  ==  'Loclinear_' ), 'method_kernel' ]  =  'Loclinear'
RES[ which( RES$method_kernel  ==  'Rejection_' ), 'method_kernel' ]  =  'Rejection'



# Guassian model:
RES_Gaus  =  RES[ which(RES$model_name == 'Gaussian'), ]
RES_Gaus   =  RES_Gaus[ , 2:7]
RES_Gaus   =  RES_Gaus[ , c(1,2,3,6)]

# Linear model:
RES_Lin   =  RES[ which(RES$model_name == 'Linear'), ]
RES_Lin   =  RES_Lin[ , 2:7]
RES_Lin   =  RES_Lin[ , c(1,2,3,6)]


RS  =  RES_Gaus  #  RES_Lin
x          =  as.character( unique( RS$dimension ) )
y          =  as.character( unique( RS$stochastic_term ) )
HT_Gaus    =  expand.grid(X=x, Y=y)
HT_Gaus$Z  =  ''

for( d in x ){
    for( st in y ){
        DF  =  RS[ intersect( which( RS$dimension        ==  d ), 
                   which( RS$stochastic_term  ==  st ) ), ]
        w = which.min( DF$MSE )
        HT_Gaus[ intersect( which( HT_Gaus$X ==  d ), 
                            which( HT_Gaus$Y  ==  st ) ), 'Z' ]  =  RS[ w, 'method_kernel']
    }
}


# colors <- c( 'green', 'yellow', "red", 'black')
# colors  =  colorRampPalette(c("red","orange","blue"),method="linear")

colors  =  colorRampPalette(c("red", "green", "blue", 'black'))(n = length( unique( HT_Gaus$Z ) ))

ggplot(HT_Gaus, aes(X, Y, fill= factor(Z) ) ) + 
    geom_tile() + 
    scale_fill_manual(values=colors) +
    xlab( 'Dimension' )  + 
    ylab( 'Stochastic term' ) + 
    theme(
        plot.title   = element_text(color="black", size=24, face="bold.italic" ),
        axis.title.x = element_text(color="black", size=24, face="bold" ),
        axis.title.y = element_text(color="black", size=24, face="bold" ),
        axis.text.x  = element_text( color="black", size=14 ),
        axis.text.y  = element_text( color="black", size=14 ),
        legend.title = element_blank(),
        legend.text = element_text( size=16, colour = "black", family = "Helvetica")
    )



--- 




# Dummy data
x <- LETTERS[1:10]
y <- paste0("var", seq(1,8))
data <- expand.grid(X=x, Y=y)
data$Z <- sample(x = c('A', 'B', 'C', 'D' ), size = 80, replace = TRUE )  # runif(40, 1, 3)



colors <- c( 'green', 'yellow', "red", 'black')

ggplot(data, aes(X, Y, fill= factor(Z) ) ) + 
    geom_tile() + 
    scale_fill_manual(values=colors)



# For continuous values of Z:
library(hrbrthemes)
ggplot(data, aes(X, Y, fill= Z)) + 
    geom_tile() +
    scale_fill_distiller(palette = "RdPu") +
    theme_ipsum()
