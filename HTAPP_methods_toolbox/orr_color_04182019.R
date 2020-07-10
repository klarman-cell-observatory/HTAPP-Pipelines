# This is Adam Haber's set of color functions.
library(colorRamps)
library(RColorBrewer)
# standardise colors

#info("Loading default colors")
default.cols = function(n){
    info(sprintf("Getting %s default colors", n))
    if(n<=20){
        #print(n)
        #info("Using 'Kelly' cols")
        kelly.cols(n)
    }else{
        warn("More than 20 requested, using 'Distinct' cols")
        distinct.cols(n)
    }
} 

wyrb.heat = colorRampPalette(c("white", "yellow3", "red2", "black"))(20)

tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


# Qualitative color schemes by Paul Tol
 tol1qualitative=c("#4477AA")
 tol2qualitative=c("#4477AA", "#CC6677")
 tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
 tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
 tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
 tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
 tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
 tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

# MONOCHROME PALETTES
# sort(brewer.pal(8,"Greens"))
redmono = c("#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0")
greenmono = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#E5F5E0", "#F7FCF5")
bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
grey8mono = c("#000000","#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
grey6mono = c("#242424", "#494949", "#6D6D6D", "#929292", "#B6B6B6", "#DBDBDB")

# EQUAL WEIGHT
# Generated with rainbow(12, s = 0.6, v = 0.75)
rainbow12equal = c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D", "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF", "#BF4DBF", "#BF4D86")
rainbow10equal = c("#BF4D4D", "#BF914D", "#A8BF4D", "#63BF4D", "#4DBF7A", "#4DBFBF", "#4D7ABF", "#634DBF", "#A84DBF", "#BF4D91")
rainbow8equal = c("#BF4D4D", "#BFA34D", "#86BF4D", "#4DBF69", "#4DBFBF", "#4D69BF", "#864DBF", "#BF4DA3")
rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
 
# Generated with package "gplots" function rich.colors(12)
rich12equal = c("#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466", "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
rich10equal = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
rich8equal = c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300")
rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
 
# Generated with package "fields" function tim.colors(12), which is said to emulate the default matlab colorset
tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
tim10equal = c("#00008F", "#0000FF", "#0070FF", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#EF0000", "#800000")
tim8equal = c("#00008F", "#0020FF", "#00AFFF", "#40FFBF", "#CFFF30", "#FF9F00", "#FF1000", "#800000")
tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
 
# Generated with sort(brewer.pal(8,"Dark2")) #Dark2, Set2
dark8equal = c("#1B9E77", "#666666", "#66A61E", "#7570B3", "#A6761D", "#D95F02", "#E6AB02", "#E7298A")
dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
set8equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#B3B3B3", "#E5C494", "#E78AC3", "#FC8D62", "#FFD92F")
set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")

## from Sam Rs compoHeatMap.R
## Get good colors for use in the heatmap. A wrapper function for
## color.palette(). Allows steps to be rescaled so that middle color
## corresponds to a given value in the range.
### ARGS:
## steps: vector of colors.
### n.steps.between: integer vector of #steps between each color given
## range.val: numeric vector of length 2 giving lower and upper limits
## of range of values that will be plotted;
### mid.val: numeric value in range.val that should be represented by
### the middle index (rounding up) in the steps vector; ignored if
### range.val is NULL.
## n.steps.final: the number of colors desired in the output vector.
### ...: Unspecified arguments are sent to color.palette().
## RETURNS:
### a vector of length n.steps.final.
get.hmap.col <- function(steps=c("blue", "cyan", "yellow", "red"), n.steps.between=c(9,1,10),
                         range.val=NULL, mid.val=NULL, n.steps.final=30,...) {
    if (!is.null(range.val)) {
        if (is.null(mid.val)) {
            mid.val=(range.val[2]-range.val[1])/2.0
        }
        mid.index=ceiling(length(n.steps.between)/2)
        ## fractional steps in the low vs. high ranges
        frac.steps.low=n.steps.between[1:mid.index]/sum(n.steps.between[1:mid.index])
        frac.steps.high=n.steps.between[(mid.index+1):length(n.steps.between)]/sum(n.steps.between[(mid.index+1):length(n.steps.between)])
        ## fraction of actual values in the low vs. high ranges
        frac.low=(mid.val-range.val[1])/(range.val[2]-range.val[1])
        frac.high=(range.val[2]-mid.val)/(range.val[2]-range.val[1])
        ## Get the right resolution and scale:
        ## n.steps is the total number of steps that will be used
        n.steps=max(255, ceiling(10^abs(log10(min(frac.low*frac.steps.low)))), ceiling(10^abs(log10(min(frac.high*frac.steps.high)))))
        ## sum(frac.high*frac.steps.high)+sum(frac.low*frac.steps.low) == 1
        n.steps.low=round(frac.low*frac.steps.low*n.steps)
        n.steps.high=round(frac.high*frac.steps.high*n.steps)
        n.steps.between=c(n.steps.low, n.steps.high)
    }
    hmcol=color.palette(steps=steps, n.steps.between=n.steps.between, ...)(n.steps.final)
    return(hmcol)
}

## Wrapper function for colorRampPalette based on
## http://stackoverflow.com/questions/13327326/r-image-function-in-r
## It allows for the definition of the number of intermediate colors
## between the main colors.  Using this option, one can stretch out
## colors that should predominate the palette spectrum. Additional
## arguments of colorRampPalette can also be added regarding the type
## and bias of the subsequent interpolation.
### ARGS:
## steps: integer.
### n.steps.between: NULL or integer.
## ...: Unspecified arguments sent to colorRampPalette().
### RETURNS:
## a color palette function, as returned by colorRampPalette.
### Usage:
## Compare pal.1 <- colorRampPalette(c("blue", "cyan", "yellow",
### "red"), bias=1)
## with
### pal.2 <- color.palette(c("blue", "cyan", "yellow", "red"),
### n.steps.between=c(10,1,10))
color.palette <- function(steps, n.steps.between=NULL, ...){
    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)
    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# https://www.materialui.co/colors
material.cols <- c("#f44336", #red
    "#E91E63", #pink
    "#9C27B0", #purple
    "#673AB7", #deep purple 
    "#3F51B5",  # indigo
    "#2196F3", # blue
    "#03A9F4", # light blue
    "#00BCD4", #cyan
    "#009688", # teal
    "#4CAF50", #green
    "#8BC34A", #light green
    "#CDDC39", # lime
    "#FFEB3B", #yellow
    "#FFC107", # amber
    "#FF9800", # organe
    "#FF5722", # deep orange
    "#795548", #brown
    "#9E9E9E", # grey
    "#607D8B" #blue grey
    )

isc.subset.cols = rev(brewer.pal(3, "Set1")) #colorRampPalette(material.700[9:12])(3)

# reverse engineered from the darksky weather app
darksky <- function(n)
{
    colorRampPalette(c(rgb(110/255, 41/255, 132/255, 1), # magenta
    rgb(33/255, 46/255, 115/255), # navy
    rgb(25/255, 96/255, 155/255), # blue
    #rgb(50/255, 150/255, 86/255), # blue 2
    rgb(80/255, 170/255, 183/255), # seafoam
    rgb(127/255, 200/255, 178/255), # teal
    rgb(235/255, 238/255, 207/255), # goldenrod
    rgb(246/255, 226/255, 155/255), # light orange
    rgb(247/255, 170/255, 86/255), # orange
    rgb(240/255, 92/255, 38/255), #scarlet
    rgb(142/255, 40/255, 11/255), #maroon
    rgb(99/255, 27/255, 7/255)) # deep red
    )(n)
}

material.heat <- function(n)
{

    mh = c(
        #"#607D8B", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        #"#03A9F4", # light blue
        "#00BCD4", #cyan
        #"#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

material.heat.new <- function(n)
{

    mh = c(
        "black",
        #"grey10",
        "grey20",
        "#2D3B79", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        "#00BCD4", #cyan
        "#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

material.700 = c("#d32f2f",
                    "#C2185B",
                    "#7B1FA2",
                    "#512DA8",
                    "#303F9F",
                    "#1976D2",
                    "#0288D1",
                    "#0097A7",
                    "#00796B",
                    "#388E3C",
                    "#689F38",
                    "#AFB42B",
                    "#FBC02D",
                    "#FFA000",
                    "#F57C00",
                    "#E64A19",
                    "#5D4037",
                    "#616161",
                    "#455A64")

material.800.heat <- function(n)
{
    m8h = c("#37474F",
    "#283593",
    "#1565C0",
    "#0277BD",
    "#00838F",
    "#00695C",
    "#2E7D32",
    "#558B2F",
    "#9E9D24",
    "#F9A825",
    "#FF8F00",
    "#EF6C00",
    "#D84315")
    colorRampPalette(m8h)(n)
}


flat.cols <- function(n)
{
    fc = c("#34495e", #wet asphalt
    "#9b59b6",  #amythest 
    "#3498db",  #peter river
    "#2ecc71",  # emerald
    #"#1abc9c",  #turquiose
    "#f1c40f",  # sunflower 
    "#e67e22",   # carrot
    "#e74c3c")   # alizarin
    #"#ecf0f1",   #clouds
    #"#95a5a6")   # concrete
    return(colorRampPalette(fc)(n))

}



heat.cols.adam <- function(n)
{
    rev(colorRampPalette(c("honeydew2", "lightgoldenrod1", "darkgoldenrod1", "firebrick1"))(n))
}

hmap.cols <- function(min.val, mid.val, max.val)
{
    get.hmap.col(mid.val=mid.val, range.val=c(min.val, max.val), 
                                steps=c("midnightblue", "blue", "cyan", "yellow", "red"), 
                                n.steps.between=c(30, 60, 60, 60))
}

## a combination of set1 and set2 with similar pinks and ugly yellow replaced 
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

brewer20 = c(brewer16, brewer.pal(12, "Set3")[c(3, 9, 8)], "violetred4")
brewer20[3] = brewer.pal(9, "Greens")[7]
brewer20[14] = brewer.pal(9, "Greens")[4]


distinct.cols <- function(n)
{
    return (intense.100 [1:n])
}

kelly.cols <- function(n)
{
    if(n <= 20)
    {
        return(kelly[1:n])
    }else
    {
        warn("Only 20 kelly colours available")
        return(kelly)
    }
}

# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# color.blind.friendly.cols <- function(n)
# {
#      if (n < length(cbbPalette))
#     {
#         return(cbbPalette[1:n])
#     }else
#     {
#         warn(sprintf("Dont have that many colourblind friendly colours. Returning %s maximally contrasting colors.", n))
#         return (c(cbbPalette[1:length(cbbPalette)], cols[1:n-length(cbbPalette)])
#     }
# }

intense.cols <-function(n)
{
    if (n < length(intense))
    {
        return(intense[1:n])
    }else if( n< length(intense.100))
    {
        return(intense.100[1:n])
    }else
    {
        warn(sprintf("Dont have that many intense colours. Returning %s maximally contrasting colors.", n))
        return (distinct.cols(n))
    }
    
}

intense.100 = c( '#6CC839','#0DA1E8','#F243F5','#3D3D12','#512257','#2BAF98','#E17516',
    '#F6A1D7','#4751C5','#D48B79','#025466','#A1AF65','#70091B','#A82F9A','#2F6401',
    '#FA3707','#A18D06','#95B0C8','#155D44','#E3A359','#F282F0','#5FC981','#C3051E',
    '#5C2D37','#993700','#E11F90','#B3C424','#94A0F5','#07326B','#0499B3','#A977F9',
    '#F6774F','#FD7C7E','#A3691D','#B33FCE','#2F3442','#C5AADA','#5485F7','#7B1E52',
    '#DC9EB2','#819722','#2176C8','#0E9164','#338D2E','#09602A','#EC61AB','#0B3D37',
    '#3B3280','#725515','#7FBAA8','#5C9A18','#611D6E','#F7A77C','#8B0938','#FB1A2A',
    '#F59A9E','#95150E','#EF5447','#35B4E5','#FC2C4F','#1D4260','#EB93E2','#6BC968',
    '#66CCA2','#1F7684','#C82EB2','#F1905E','#ED127F','#2B4A2F','#7EBCC0','#644037',
    '#ACB4DE','#81C1E0','#CCA660','#5B1E3C','#FA3774','#81A200','#7C79EB','#87C875',
    '#8F1121','#473142','#DA5BEE','#751D61','#EC63D9','#437C04','#9059DF','#A3B235',
    '#56450F','#13A367','#E47EFC','#E1AD4F','#236E42','#3875B8','#1B5118','#3F490F',
    '#F1153C','#5A2F29','#0F6050','#1E6976')

# http://tools.medialab.sciences-po.fr/iwanthue/
intense = c("#F41AA7",
    "#12802C",
    "#F56700",
    "#4371B9",
    "#5C270E",
    "#B9BD09",
    "#D99B68",
    "#157D71",
    "#40334F",
    "#D80F35",
    "#CFA6DA",
    "#B8085C",
    "#B633C3",
    "#3D4326",
    "#6447AD",
    "#76B57E",
    "#F5A83E",
    "#F570BE",
    "#851122",
    "#31360C")

cols = c(
        "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
        
        "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
        "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
        "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
        "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
        "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
        "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
        "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
        "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",

        "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
        "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
        "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
        "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
        "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
        "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
        "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
        "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
        
        "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
        "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B")

kelly = c(
    "#00538A", # Strong Blue
    "#C10020", # Vivid Red
    "#007D34", # Vivid Green
    "#FFB300", # Vivid Yellow
    "#803E75", # Strong Purple
    "#FF6800", # Vivid Orange
    "#A6BDD7", # Very Light Blue
    "#CEA262", # Grayish Yellow
    "#817066", # Medium Gray
    "#F6768E", # Strong Purplish Pink
    "#FF7A5C", # Strong Yellowish Pink
    "#53377A", # Strong Violet
    "#FF8E00", # Vivid Orange Yellow
    "#B32851", # Strong Purplish Red
    "#F4C800", # Vivid Greenish Yellow
    "#7F180D", # Strong Reddish Brown
    "#93AA00", # Vivid Yellowish Green
    "#593315", # Deep Yellowish Brown
    "#F13A13", # Vivid Reddish Orange
    "#232C16" # Dark Olive Green
    )

color.kelly <- function(n)
{
   if(n < length(kelly))
   {
	return (sample(kelly, n))
   }else{
	return (colorRampPalette(kelly)(n))
   }
}

# plot some colors
 pal <- function(col, border = "light gray", ...)
 {
     n <- length(col)
     plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "", ...)
     rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
 }

 rainbow.hcl <- function(n=10)
 {
    library(colorspace)
    rainbow_hcl(n, c = 50, l = 70, start = 0, end = 360*(n-1)/n,
        gamma = NULL, fixup = TRUE, alpha = 1)
 }

 heat.hcl <- function(n=10)
 {
    library(colorspace)
    rev(heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 2)))
 }

intense.16 = c("#FFBC6B", "#727EFF", "#8FE300", "#FF075F", "#01EB50", "#9A5D9E", "#F1DE3B", "#85B6FF", "#9A6A0F", "#25D0FF", "#BB4F67", "#00E59F", "#FFC6C7", "#00873A", "#00A1B2",
"#2D8163")
 

intense.16.2 = toupper(c("#bfca41", "#ab48f9", "#019f25", "#b5318e", "#bf9600", "#0190f3", "#ff8e39", "#eb98ff", "#007a54", "#ff649e", "#94ccc6", "#ca232d", "#bcc0da", "#a9466f", "#bfc68f", "#81624d"))

ad.cubehelix.old <- function(n)
{
    library(rje)
    cubeHelix(n, gamma = 0.4, hue=1.6, r=-0.9, start = -3)[ceiling(n/10):n]
}

ad.cubehelix <- function(n)
{
    library(rje)
    cubeHelix(n, gamma = 0.4, hue=1.6, r=-0.9, start = -3)[ceiling(n/10):floor(n*9/10)]
}

distinct.phrogz <- function(n)
{
    x = c("#BF3030","#FFA280","#A68A53","#EEF2B6","#00BF80","#206C80","#262A33","#7340FF","#D9A3CE","#331A1A","#B2622D","#FFF240","#829926","#2D594A","#0088CC",
    "#535EA6","#655673","#E63995","#A67C7C","#593E2D","#333226","#41F200","#7CA698","#002999","#BFC8FF","#B630BF","#592D3E","#F24100","#F29D3D","#555916",
    "#24661A","#3DE6F2","#101D40","#110040","#6D1D73","#66001B")
    return(x[1:n])
}

# from python package 'palettable'
cubehelix.classic16 = c('#000000', '#160A22', '#182044', '#103E53', 
    '#0E5E4A', '#237433', '#507D23', '#8A7A2D', '#BE7555', '#DA7991', 
    '#DB8ACB', '#CCA7F0', '#BFC9FB', '#C3E5F4', '#DCF6EF', '#FFFFFF')

cubehelix1.16 = c('#000000', '#1B0F00', '#411704', '#681B20', 
    '#85214B', '#932D7E', '#9042AF', '#8160D2', '#6F83E3', 
    '#63A6E2', '#65C5D3', '#78DBC2', '#99E9B9', '#C1F0BF', '#E6F5D8', '#FFFFFF')

cubehelix2.16 = c('#000000', '#001C0E', '#00332F', '#07415B', 
    '#234787', '#4E48A8', '#8148B8', '#B14DB5', '#D65AA5', 
    '#EB718F', '#EE8E80', '#E6AF7F', '#DBCE90', '#D8E7B2', '#E2F7DB', '#FFFFFF')

jim_special_16 = c('#000000', '#160A22', '#251944', '#2F2B63', 
    '#34417D', '#375892', '#3B70A0', '#4089A9', '#4AA0AD', '#59B5AF',
     '#6DC7B1', '#86D6B4', '#A3E3BD', '#C3EDCB', '#E2F6E1', '#FFFFFF')

perceptual_rainbow_12 = c('#873B61', '#8F489D',
                          '#7966CF', '#677BDC', '#5492DF', '#45AAD7', '#3BC0C5', '#3CD2AC',
                          '#47DF91', '#5DE578', '#A1E35F', '#E9D575')

perceptual_rainbow_16 = c('#873B61', '#8F407F', '#8F489D', '#8755B9', 
    '#7966CF', '#677BDC', '#5492DF', '#45AAD7', '#3BC0C5', '#3CD2AC', 
    '#47DF91', '#5DE578', '#7CE767', '#A1E35F', '#C6DC64', '#E9D575')

red_16 = c('#000000', '#130C23', '#2C1641', '#49205A', 
    '#68296B', '#863576', '#A2437C', '#B9537E', '#CC6680', 
    '#D87C82', '#E19488', '#E5AC93', '#E8C4A4', '#ECDBBD', '#F2EEDB', '#FFFFFF')


matplotlib.magma = c(rgb(0.001462, 0.000466, 0.013866),
               rgb(0.002258, 0.001295, 0.018331),
               rgb(0.003279, 0.002305, 0.023708),
               rgb(0.004512, 0.003490, 0.029965),
               rgb(0.005950, 0.004843, 0.037130),
               rgb(0.007588, 0.006356, 0.044973),
               rgb(0.009426, 0.008022, 0.052844),
               rgb(0.011465, 0.009828, 0.060750),
               rgb(0.013708, 0.011771, 0.068667),
               rgb(0.016156, 0.013840, 0.076603),
               rgb(0.018815, 0.016026, 0.084584),
               rgb(0.021692, 0.018320, 0.092610),
               rgb(0.024792, 0.020715, 0.100676),
               rgb(0.028123, 0.023201, 0.108787),
               rgb(0.031696, 0.025765, 0.116965),
               rgb(0.035520, 0.028397, 0.125209),
               rgb(0.039608, 0.031090, 0.133515),
               rgb(0.043830, 0.033830, 0.141886),
               rgb(0.048062, 0.036607, 0.150327),
               rgb(0.052320, 0.039407, 0.158841),
               rgb(0.056615, 0.042160, 0.167446),
               rgb(0.060949, 0.044794, 0.176129),
               rgb(0.065330, 0.047318, 0.184892),
               rgb(0.069764, 0.049726, 0.193735),
               rgb(0.074257, 0.052017, 0.202660),
               rgb(0.078815, 0.054184, 0.211667),
               rgb(0.083446, 0.056225, 0.220755),
               rgb(0.088155, 0.058133, 0.229922),
               rgb(0.092949, 0.059904, 0.239164),
               rgb(0.097833, 0.061531, 0.248477),
               rgb(0.102815, 0.063010, 0.257854),
               rgb(0.107899, 0.064335, 0.267289),
               rgb(0.113094, 0.065492, 0.276784),
               rgb(0.118405, 0.066479, 0.286321),
               rgb(0.123833, 0.067295, 0.295879),
               rgb(0.129380, 0.067935, 0.305443),
               rgb(0.135053, 0.068391, 0.315000),
               rgb(0.140858, 0.068654, 0.324538),
               rgb(0.146785, 0.068738, 0.334011),
               rgb(0.152839, 0.068637, 0.343404),
               rgb(0.159018, 0.068354, 0.352688),
               rgb(0.165308, 0.067911, 0.361816),
               rgb(0.171713, 0.067305, 0.370771),
               rgb(0.178212, 0.066576, 0.379497),
               rgb(0.184801, 0.065732, 0.387973),
               rgb(0.191460, 0.064818, 0.396152),
               rgb(0.198177, 0.063862, 0.404009),
               rgb(0.204935, 0.062907, 0.411514),
               rgb(0.211718, 0.061992, 0.418647),
               rgb(0.218512, 0.061158, 0.425392),
               rgb(0.225302, 0.060445, 0.431742),
               rgb(0.232077, 0.059889, 0.437695),
               rgb(0.238826, 0.059517, 0.443256),
               rgb(0.245543, 0.059352, 0.448436),
               rgb(0.252220, 0.059415, 0.453248),
               rgb(0.258857, 0.059706, 0.457710),
               rgb(0.265447, 0.060237, 0.461840),
               rgb(0.271994, 0.060994, 0.465660),
               rgb(0.278493, 0.061978, 0.469190),
               rgb(0.284951, 0.063168, 0.472451),
               rgb(0.291366, 0.064553, 0.475462),
               rgb(0.297740, 0.066117, 0.478243),
               rgb(0.304081, 0.067835, 0.480812),
               rgb(0.310382, 0.069702, 0.483186),
               rgb(0.316654, 0.071690, 0.485380),
               rgb(0.322899, 0.073782, 0.487408),
               rgb(0.329114, 0.075972, 0.489287),
               rgb(0.335308, 0.078236, 0.491024),
               rgb(0.341482, 0.080564, 0.492631),
               rgb(0.347636, 0.082946, 0.494121),
               rgb(0.353773, 0.085373, 0.495501),
               rgb(0.359898, 0.087831, 0.496778),
               rgb(0.366012, 0.090314, 0.497960),
               rgb(0.372116, 0.092816, 0.499053),
               rgb(0.378211, 0.095332, 0.500067),
               rgb(0.384299, 0.097855, 0.501002),
               rgb(0.390384, 0.100379, 0.501864),
               rgb(0.396467, 0.102902, 0.502658),
               rgb(0.402548, 0.105420, 0.503386),
               rgb(0.408629, 0.107930, 0.504052),
               rgb(0.414709, 0.110431, 0.504662),
               rgb(0.420791, 0.112920, 0.505215),
               rgb(0.426877, 0.115395, 0.505714),
               rgb(0.432967, 0.117855, 0.506160),
               rgb(0.439062, 0.120298, 0.506555),
               rgb(0.445163, 0.122724, 0.506901),
               rgb(0.451271, 0.125132, 0.507198),
               rgb(0.457386, 0.127522, 0.507448),
               rgb(0.463508, 0.129893, 0.507652),
               rgb(0.469640, 0.132245, 0.507809),
               rgb(0.475780, 0.134577, 0.507921),
               rgb(0.481929, 0.136891, 0.507989),
               rgb(0.488088, 0.139186, 0.508011),
               rgb(0.494258, 0.141462, 0.507988),
               rgb(0.500438, 0.143719, 0.507920),
               rgb(0.506629, 0.145958, 0.507806),
               rgb(0.512831, 0.148179, 0.507648),
               rgb(0.519045, 0.150383, 0.507443),
               rgb(0.525270, 0.152569, 0.507192),
               rgb(0.531507, 0.154739, 0.506895),
               rgb(0.537755, 0.156894, 0.506551),
               rgb(0.544015, 0.159033, 0.506159),
               rgb(0.550287, 0.161158, 0.505719),
               rgb(0.556571, 0.163269, 0.505230),
               rgb(0.562866, 0.165368, 0.504692),
               rgb(0.569172, 0.167454, 0.504105),
               rgb(0.575490, 0.169530, 0.503466),
               rgb(0.581819, 0.171596, 0.502777),
               rgb(0.588158, 0.173652, 0.502035),
               rgb(0.594508, 0.175701, 0.501241),
               rgb(0.600868, 0.177743, 0.500394),
               rgb(0.607238, 0.179779, 0.499492),
               rgb(0.613617, 0.181811, 0.498536),
               rgb(0.620005, 0.183840, 0.497524),
               rgb(0.626401, 0.185867, 0.496456),
               rgb(0.632805, 0.187893, 0.495332),
               rgb(0.639216, 0.189921, 0.494150),
               rgb(0.645633, 0.191952, 0.492910),
               rgb(0.652056, 0.193986, 0.491611),
               rgb(0.658483, 0.196027, 0.490253),
               rgb(0.664915, 0.198075, 0.488836),
               rgb(0.671349, 0.200133, 0.487358),
               rgb(0.677786, 0.202203, 0.485819),
               rgb(0.684224, 0.204286, 0.484219),
               rgb(0.690661, 0.206384, 0.482558),
               rgb(0.697098, 0.208501, 0.480835),
               rgb(0.703532, 0.210638, 0.479049),
               rgb(0.709962, 0.212797, 0.477201),
               rgb(0.716387, 0.214982, 0.475290),
               rgb(0.722805, 0.217194, 0.473316),
               rgb(0.729216, 0.219437, 0.471279),
               rgb(0.735616, 0.221713, 0.469180),
               rgb(0.742004, 0.224025, 0.467018),
               rgb(0.748378, 0.226377, 0.464794),
               rgb(0.754737, 0.228772, 0.462509),
               rgb(0.761077, 0.231214, 0.460162),
               rgb(0.767398, 0.233705, 0.457755),
               rgb(0.773695, 0.236249, 0.455289),
               rgb(0.779968, 0.238851, 0.452765),
               rgb(0.786212, 0.241514, 0.450184),
               rgb(0.792427, 0.244242, 0.447543),
               rgb(0.798608, 0.247040, 0.444848),
               rgb(0.804752, 0.249911, 0.442102),
               rgb(0.810855, 0.252861, 0.439305),
               rgb(0.816914, 0.255895, 0.436461),
               rgb(0.822926, 0.259016, 0.433573),
               rgb(0.828886, 0.262229, 0.430644),
               rgb(0.834791, 0.265540, 0.427671),
               rgb(0.840636, 0.268953, 0.424666),
               rgb(0.846416, 0.272473, 0.421631),
               rgb(0.852126, 0.276106, 0.418573),
               rgb(0.857763, 0.279857, 0.415496),
               rgb(0.863320, 0.283729, 0.412403),
               rgb(0.868793, 0.287728, 0.409303),
               rgb(0.874176, 0.291859, 0.406205),
               rgb(0.879464, 0.296125, 0.403118),
               rgb(0.884651, 0.300530, 0.400047),
               rgb(0.889731, 0.305079, 0.397002),
               rgb(0.894700, 0.309773, 0.393995),
               rgb(0.899552, 0.314616, 0.391037),
               rgb(0.904281, 0.319610, 0.388137),
               rgb(0.908884, 0.324755, 0.385308),
               rgb(0.913354, 0.330052, 0.382563),
               rgb(0.917689, 0.335500, 0.379915),
               rgb(0.921884, 0.341098, 0.377376),
               rgb(0.925937, 0.346844, 0.374959),
               rgb(0.929845, 0.352734, 0.372677),
               rgb(0.933606, 0.358764, 0.370541),
               rgb(0.937221, 0.364929, 0.368567),
               rgb(0.940687, 0.371224, 0.366762),
               rgb(0.944006, 0.377643, 0.365136),
               rgb(0.947180, 0.384178, 0.363701),
               rgb(0.950210, 0.390820, 0.362468),
               rgb(0.953099, 0.397563, 0.361438),
               rgb(0.955849, 0.404400, 0.360619),
               rgb(0.958464, 0.411324, 0.360014),
               rgb(0.960949, 0.418323, 0.359630),
               rgb(0.963310, 0.425390, 0.359469),
               rgb(0.965549, 0.432519, 0.359529),
               rgb(0.967671, 0.439703, 0.359810),
               rgb(0.969680, 0.446936, 0.360311),
               rgb(0.971582, 0.454210, 0.361030),
               rgb(0.973381, 0.461520, 0.361965),
               rgb(0.975082, 0.468861, 0.363111),
               rgb(0.976690, 0.476226, 0.364466),
               rgb(0.978210, 0.483612, 0.366025),
               rgb(0.979645, 0.491014, 0.367783),
               rgb(0.981000, 0.498428, 0.369734),
               rgb(0.982279, 0.505851, 0.371874),
               rgb(0.983485, 0.513280, 0.374198),
               rgb(0.984622, 0.520713, 0.376698),
               rgb(0.985693, 0.528148, 0.379371),
               rgb(0.986700, 0.535582, 0.382210),
               rgb(0.987646, 0.543015, 0.385210),
               rgb(0.988533, 0.550446, 0.388365),
               rgb(0.989363, 0.557873, 0.391671),
               rgb(0.990138, 0.565296, 0.395122),
               rgb(0.990871, 0.572706, 0.398714),
               rgb(0.991558, 0.580107, 0.402441),
               rgb(0.992196, 0.587502, 0.406299),
               rgb(0.992785, 0.594891, 0.410283),
               rgb(0.993326, 0.602275, 0.414390),
               rgb(0.993834, 0.609644, 0.418613),
               rgb(0.994309, 0.616999, 0.422950),
               rgb(0.994738, 0.624350, 0.427397),
               rgb(0.995122, 0.631696, 0.431951),
               rgb(0.995480, 0.639027, 0.436607),
               rgb(0.995810, 0.646344, 0.441361),
               rgb(0.996096, 0.653659, 0.446213),
               rgb(0.996341, 0.660969, 0.451160),
               rgb(0.996580, 0.668256, 0.456192),
               rgb(0.996775, 0.675541, 0.461314),
               rgb(0.996925, 0.682828, 0.466526),
               rgb(0.997077, 0.690088, 0.471811),
               rgb(0.997186, 0.697349, 0.477182),
               rgb(0.997254, 0.704611, 0.482635),
               rgb(0.997325, 0.711848, 0.488154),
               rgb(0.997351, 0.719089, 0.493755),
               rgb(0.997351, 0.726324, 0.499428),
               rgb(0.997341, 0.733545, 0.505167),
               rgb(0.997285, 0.740772, 0.510983),
               rgb(0.997228, 0.747981, 0.516859),
               rgb(0.997138, 0.755190, 0.522806),
               rgb(0.997019, 0.762398, 0.528821),
               rgb(0.996898, 0.769591, 0.534892),
               rgb(0.996727, 0.776795, 0.541039),
               rgb(0.996571, 0.783977, 0.547233),
               rgb(0.996369, 0.791167, 0.553499),
               rgb(0.996162, 0.798348, 0.559820),
               rgb(0.995932, 0.805527, 0.566202),
               rgb(0.995680, 0.812706, 0.572645),
               rgb(0.995424, 0.819875, 0.579140),
               rgb(0.995131, 0.827052, 0.585701),
               rgb(0.994851, 0.834213, 0.592307),
               rgb(0.994524, 0.841387, 0.598983),
               rgb(0.994222, 0.848540, 0.605696),
               rgb(0.993866, 0.855711, 0.612482),
               rgb(0.993545, 0.862859, 0.619299),
               rgb(0.993170, 0.870024, 0.626189),
               rgb(0.992831, 0.877168, 0.633109),
               rgb(0.992440, 0.884330, 0.640099),
               rgb(0.992089, 0.891470, 0.647116),
               rgb(0.991688, 0.898627, 0.654202),
               rgb(0.991332, 0.905763, 0.661309),
               rgb(0.990930, 0.912915, 0.668481),
               rgb(0.990570, 0.920049, 0.675675),
               rgb(0.990175, 0.927196, 0.682926),
               rgb(0.989815, 0.934329, 0.690198),
               rgb(0.989434, 0.941470, 0.697519),
               rgb(0.989077, 0.948604, 0.704863),
               rgb(0.988717, 0.955742, 0.712242),
               rgb(0.988367, 0.962878, 0.719649),
               rgb(0.988033, 0.970012, 0.727077),
               rgb(0.987691, 0.977154, 0.734536),
               rgb(0.987387, 0.984288, 0.742002),
               rgb(0.987053, 0.991438, 0.749504))


matplotlib.inferno <- c(rgb(0.001462, 0.000466, 0.013866),
                 rgb(0.002267, 0.001270, 0.018570),
                 rgb(0.003299, 0.002249, 0.024239),
                 rgb(0.004547, 0.003392, 0.030909),
                 rgb(0.006006, 0.004692, 0.038558),
                 rgb(0.007676, 0.006136, 0.046836),
                 rgb(0.009561, 0.007713, 0.055143),
                 rgb(0.011663, 0.009417, 0.063460),
                 rgb(0.013995, 0.011225, 0.071862),
                 rgb(0.016561, 0.013136, 0.080282),
                 rgb(0.019373, 0.015133, 0.088767),
                 rgb(0.022447, 0.017199, 0.097327),
                 rgb(0.025793, 0.019331, 0.105930),
                 rgb(0.029432, 0.021503, 0.114621),
                 rgb(0.033385, 0.023702, 0.123397),
                 rgb(0.037668, 0.025921, 0.132232),
                 rgb(0.042253, 0.028139, 0.141141),
                 rgb(0.046915, 0.030324, 0.150164),
                 rgb(0.051644, 0.032474, 0.159254),
                 rgb(0.056449, 0.034569, 0.168414),
                 rgb(0.061340, 0.036590, 0.177642),
                 rgb(0.066331, 0.038504, 0.186962),
                 rgb(0.071429, 0.040294, 0.196354),
                 rgb(0.076637, 0.041905, 0.205799),
                 rgb(0.081962, 0.043328, 0.215289),
                 rgb(0.087411, 0.044556, 0.224813),
                 rgb(0.092990, 0.045583, 0.234358),
                 rgb(0.098702, 0.046402, 0.243904),
                 rgb(0.104551, 0.047008, 0.253430),
                 rgb(0.110536, 0.047399, 0.262912),
                 rgb(0.116656, 0.047574, 0.272321),
                 rgb(0.122908, 0.047536, 0.281624),
                 rgb(0.129285, 0.047293, 0.290788),
                 rgb(0.135778, 0.046856, 0.299776),
                 rgb(0.142378, 0.046242, 0.308553),
                 rgb(0.149073, 0.045468, 0.317085),
                 rgb(0.155850, 0.044559, 0.325338),
                 rgb(0.162689, 0.043554, 0.333277),
                 rgb(0.169575, 0.042489, 0.340874),
                 rgb(0.176493, 0.041402, 0.348111),
                 rgb(0.183429, 0.040329, 0.354971),
                 rgb(0.190367, 0.039309, 0.361447),
                 rgb(0.197297, 0.038400, 0.367535),
                 rgb(0.204209, 0.037632, 0.373238),
                 rgb(0.211095, 0.037030, 0.378563),
                 rgb(0.217949, 0.036615, 0.383522),
                 rgb(0.224763, 0.036405, 0.388129),
                 rgb(0.231538, 0.036405, 0.392400),
                 rgb(0.238273, 0.036621, 0.396353),
                 rgb(0.244967, 0.037055, 0.400007),
                 rgb(0.251620, 0.037705, 0.403378),
                 rgb(0.258234, 0.038571, 0.406485),
                 rgb(0.264810, 0.039647, 0.409345),
                 rgb(0.271347, 0.040922, 0.411976),
                 rgb(0.277850, 0.042353, 0.414392),
                 rgb(0.284321, 0.043933, 0.416608),
                 rgb(0.290763, 0.045644, 0.418637),
                 rgb(0.297178, 0.047470, 0.420491),
                 rgb(0.303568, 0.049396, 0.422182),
                 rgb(0.309935, 0.051407, 0.423721),
                 rgb(0.316282, 0.053490, 0.425116),
                 rgb(0.322610, 0.055634, 0.426377),
                 rgb(0.328921, 0.057827, 0.427511),
                 rgb(0.335217, 0.060060, 0.428524),
                 rgb(0.341500, 0.062325, 0.429425),
                 rgb(0.347771, 0.064616, 0.430217),
                 rgb(0.354032, 0.066925, 0.430906),
                 rgb(0.360284, 0.069247, 0.431497),
                 rgb(0.366529, 0.071579, 0.431994),
                 rgb(0.372768, 0.073915, 0.432400),
                 rgb(0.379001, 0.076253, 0.432719),
                 rgb(0.385228, 0.078591, 0.432955),
                 rgb(0.391453, 0.080927, 0.433109),
                 rgb(0.397674, 0.083257, 0.433183),
                 rgb(0.403894, 0.085580, 0.433179),
                 rgb(0.410113, 0.087896, 0.433098),
                 rgb(0.416331, 0.090203, 0.432943),
                 rgb(0.422549, 0.092501, 0.432714),
                 rgb(0.428768, 0.094790, 0.432412),
                 rgb(0.434987, 0.097069, 0.432039),
                 rgb(0.441207, 0.099338, 0.431594),
                 rgb(0.447428, 0.101597, 0.431080),
                 rgb(0.453651, 0.103848, 0.430498),
                 rgb(0.459875, 0.106089, 0.429846),
                 rgb(0.466100, 0.108322, 0.429125),
                 rgb(0.472328, 0.110547, 0.428334),
                 rgb(0.478558, 0.112764, 0.427475),
                 rgb(0.484789, 0.114974, 0.426548),
                 rgb(0.491022, 0.117179, 0.425552),
                 rgb(0.497257, 0.119379, 0.424488),
                 rgb(0.503493, 0.121575, 0.423356),
                 rgb(0.509730, 0.123769, 0.422156),
                 rgb(0.515967, 0.125960, 0.420887),
                 rgb(0.522206, 0.128150, 0.419549),
                 rgb(0.528444, 0.130341, 0.418142),
                 rgb(0.534683, 0.132534, 0.416667),
                 rgb(0.540920, 0.134729, 0.415123),
                 rgb(0.547157, 0.136929, 0.413511),
                 rgb(0.553392, 0.139134, 0.411829),
                 rgb(0.559624, 0.141346, 0.410078),
                 rgb(0.565854, 0.143567, 0.408258),
                 rgb(0.572081, 0.145797, 0.406369),
                 rgb(0.578304, 0.148039, 0.404411),
                 rgb(0.584521, 0.150294, 0.402385),
                 rgb(0.590734, 0.152563, 0.400290),
                 rgb(0.596940, 0.154848, 0.398125),
                 rgb(0.603139, 0.157151, 0.395891),
                 rgb(0.609330, 0.159474, 0.393589),
                 rgb(0.615513, 0.161817, 0.391219),
                 rgb(0.621685, 0.164184, 0.388781),
                 rgb(0.627847, 0.166575, 0.386276),
                 rgb(0.633998, 0.168992, 0.383704),
                 rgb(0.640135, 0.171438, 0.381065),
                 rgb(0.646260, 0.173914, 0.378359),
                 rgb(0.652369, 0.176421, 0.375586),
                 rgb(0.658463, 0.178962, 0.372748),
                 rgb(0.664540, 0.181539, 0.369846),
                 rgb(0.670599, 0.184153, 0.366879),
                 rgb(0.676638, 0.186807, 0.363849),
                 rgb(0.682656, 0.189501, 0.360757),
                 rgb(0.688653, 0.192239, 0.357603),
                 rgb(0.694627, 0.195021, 0.354388),
                 rgb(0.700576, 0.197851, 0.351113),
                 rgb(0.706500, 0.200728, 0.347777),
                 rgb(0.712396, 0.203656, 0.344383),
                 rgb(0.718264, 0.206636, 0.340931),
                 rgb(0.724103, 0.209670, 0.337424),
                 rgb(0.729909, 0.212759, 0.333861),
                 rgb(0.735683, 0.215906, 0.330245),
                 rgb(0.741423, 0.219112, 0.326576),
                 rgb(0.747127, 0.222378, 0.322856),
                 rgb(0.752794, 0.225706, 0.319085),
                 rgb(0.758422, 0.229097, 0.315266),
                 rgb(0.764010, 0.232554, 0.311399),
                 rgb(0.769556, 0.236077, 0.307485),
                 rgb(0.775059, 0.239667, 0.303526),
                 rgb(0.780517, 0.243327, 0.299523),
                 rgb(0.785929, 0.247056, 0.295477),
                 rgb(0.791293, 0.250856, 0.291390),
                 rgb(0.796607, 0.254728, 0.287264),
                 rgb(0.801871, 0.258674, 0.283099),
                 rgb(0.807082, 0.262692, 0.278898),
                 rgb(0.812239, 0.266786, 0.274661),
                 rgb(0.817341, 0.270954, 0.270390),
                 rgb(0.822386, 0.275197, 0.266085),
                 rgb(0.827372, 0.279517, 0.261750),
                 rgb(0.832299, 0.283913, 0.257383),
                 rgb(0.837165, 0.288385, 0.252988),
                 rgb(0.841969, 0.292933, 0.248564),
                 rgb(0.846709, 0.297559, 0.244113),
                 rgb(0.851384, 0.302260, 0.239636),
                 rgb(0.855992, 0.307038, 0.235133),
                 rgb(0.860533, 0.311892, 0.230606),
                 rgb(0.865006, 0.316822, 0.226055),
                 rgb(0.869409, 0.321827, 0.221482),
                 rgb(0.873741, 0.326906, 0.216886),
                 rgb(0.878001, 0.332060, 0.212268),
                 rgb(0.882188, 0.337287, 0.207628),
                 rgb(0.886302, 0.342586, 0.202968),
                 rgb(0.890341, 0.347957, 0.198286),
                 rgb(0.894305, 0.353399, 0.193584),
                 rgb(0.898192, 0.358911, 0.188860),
                 rgb(0.902003, 0.364492, 0.184116),
                 rgb(0.905735, 0.370140, 0.179350),
                 rgb(0.909390, 0.375856, 0.174563),
                 rgb(0.912966, 0.381636, 0.169755),
                 rgb(0.916462, 0.387481, 0.164924),
                 rgb(0.919879, 0.393389, 0.160070),
                 rgb(0.923215, 0.399359, 0.155193),
                 rgb(0.926470, 0.405389, 0.150292),
                 rgb(0.929644, 0.411479, 0.145367),
                 rgb(0.932737, 0.417627, 0.140417),
                 rgb(0.935747, 0.423831, 0.135440),
                 rgb(0.938675, 0.430091, 0.130438),
                 rgb(0.941521, 0.436405, 0.125409),
                 rgb(0.944285, 0.442772, 0.120354),
                 rgb(0.946965, 0.449191, 0.115272),
                 rgb(0.949562, 0.455660, 0.110164),
                 rgb(0.952075, 0.462178, 0.105031),
                 rgb(0.954506, 0.468744, 0.099874),
                 rgb(0.956852, 0.475356, 0.094695),
                 rgb(0.959114, 0.482014, 0.089499),
                 rgb(0.961293, 0.488716, 0.084289),
                 rgb(0.963387, 0.495462, 0.079073),
                 rgb(0.965397, 0.502249, 0.073859),
                 rgb(0.967322, 0.509078, 0.068659),
                 rgb(0.969163, 0.515946, 0.063488),
                 rgb(0.970919, 0.522853, 0.058367),
                 rgb(0.972590, 0.529798, 0.053324),
                 rgb(0.974176, 0.536780, 0.048392),
                 rgb(0.975677, 0.543798, 0.043618),
                 rgb(0.977092, 0.550850, 0.039050),
                 rgb(0.978422, 0.557937, 0.034931),
                 rgb(0.979666, 0.565057, 0.031409),
                 rgb(0.980824, 0.572209, 0.028508),
                 rgb(0.981895, 0.579392, 0.026250),
                 rgb(0.982881, 0.586606, 0.024661),
                 rgb(0.983779, 0.593849, 0.023770),
                 rgb(0.984591, 0.601122, 0.023606),
                 rgb(0.985315, 0.608422, 0.024202),
                 rgb(0.985952, 0.615750, 0.025592),
                 rgb(0.986502, 0.623105, 0.027814),
                 rgb(0.986964, 0.630485, 0.030908),
                 rgb(0.987337, 0.637890, 0.034916),
                 rgb(0.987622, 0.645320, 0.039886),
                 rgb(0.987819, 0.652773, 0.045581),
                 rgb(0.987926, 0.660250, 0.051750),
                 rgb(0.987945, 0.667748, 0.058329),
                 rgb(0.987874, 0.675267, 0.065257),
                 rgb(0.987714, 0.682807, 0.072489),
                 rgb(0.987464, 0.690366, 0.079990),
                 rgb(0.987124, 0.697944, 0.087731),
                 rgb(0.986694, 0.705540, 0.095694),
                 rgb(0.986175, 0.713153, 0.103863),
                 rgb(0.985566, 0.720782, 0.112229),
                 rgb(0.984865, 0.728427, 0.120785),
                 rgb(0.984075, 0.736087, 0.129527),
                 rgb(0.983196, 0.743758, 0.138453),
                 rgb(0.982228, 0.751442, 0.147565),
                 rgb(0.981173, 0.759135, 0.156863),
                 rgb(0.980032, 0.766837, 0.166353),
                 rgb(0.978806, 0.774545, 0.176037),
                 rgb(0.977497, 0.782258, 0.185923),
                 rgb(0.976108, 0.789974, 0.196018),
                 rgb(0.974638, 0.797692, 0.206332),
                 rgb(0.973088, 0.805409, 0.216877),
                 rgb(0.971468, 0.813122, 0.227658),
                 rgb(0.969783, 0.820825, 0.238686),
                 rgb(0.968041, 0.828515, 0.249972),
                 rgb(0.966243, 0.836191, 0.261534),
                 rgb(0.964394, 0.843848, 0.273391),
                 rgb(0.962517, 0.851476, 0.285546),
                 rgb(0.960626, 0.859069, 0.298010),
                 rgb(0.958720, 0.866624, 0.310820),
                 rgb(0.956834, 0.874129, 0.323974),
                 rgb(0.954997, 0.881569, 0.337475),
                 rgb(0.953215, 0.888942, 0.351369),
                 rgb(0.951546, 0.896226, 0.365627),
                 rgb(0.950018, 0.903409, 0.380271),
                 rgb(0.948683, 0.910473, 0.395289),
                 rgb(0.947594, 0.917399, 0.410665),
                 rgb(0.946809, 0.924168, 0.426373),
                 rgb(0.946392, 0.930761, 0.442367),
                 rgb(0.946403, 0.937159, 0.458592),
                 rgb(0.946903, 0.943348, 0.474970),
                 rgb(0.947937, 0.949318, 0.491426),
                 rgb(0.949545, 0.955063, 0.507860),
                 rgb(0.951740, 0.960587, 0.524203),
                 rgb(0.954529, 0.965896, 0.540361),
                 rgb(0.957896, 0.971003, 0.556275),
                 rgb(0.961812, 0.975924, 0.571925),
                 rgb(0.966249, 0.980678, 0.587206),
                 rgb(0.971162, 0.985282, 0.602154),
                 rgb(0.976511, 0.989753, 0.616760),
                 rgb(0.982257, 0.994109, 0.631017),
                 rgb(0.988362, 0.998364, 0.644924))

matplotlib.viridis <- c(rgb(0.267004, 0.004874, 0.329415),
                 rgb(0.268510, 0.009605, 0.335427),
                 rgb(0.269944, 0.014625, 0.341379),
                 rgb(0.271305, 0.019942, 0.347269),
                 rgb(0.272594, 0.025563, 0.353093),
                 rgb(0.273809, 0.031497, 0.358853),
                 rgb(0.274952, 0.037752, 0.364543),
                 rgb(0.276022, 0.044167, 0.370164),
                 rgb(0.277018, 0.050344, 0.375715),
                 rgb(0.277941, 0.056324, 0.381191),
                 rgb(0.278791, 0.062145, 0.386592),
                 rgb(0.279566, 0.067836, 0.391917),
                 rgb(0.280267, 0.073417, 0.397163),
                 rgb(0.280894, 0.078907, 0.402329),
                 rgb(0.281446, 0.084320, 0.407414),
                 rgb(0.281924, 0.089666, 0.412415),
                 rgb(0.282327, 0.094955, 0.417331),
                 rgb(0.282656, 0.100196, 0.422160),
                 rgb(0.282910, 0.105393, 0.426902),
                 rgb(0.283091, 0.110553, 0.431554),
                 rgb(0.283197, 0.115680, 0.436115),
                 rgb(0.283229, 0.120777, 0.440584),
                 rgb(0.283187, 0.125848, 0.444960),
                 rgb(0.283072, 0.130895, 0.449241),
                 rgb(0.282884, 0.135920, 0.453427),
                 rgb(0.282623, 0.140926, 0.457517),
                 rgb(0.282290, 0.145912, 0.461510),
                 rgb(0.281887, 0.150881, 0.465405),
                 rgb(0.281412, 0.155834, 0.469201),
                 rgb(0.280868, 0.160771, 0.472899),
                 rgb(0.280255, 0.165693, 0.476498),
                 rgb(0.279574, 0.170599, 0.479997),
                 rgb(0.278826, 0.175490, 0.483397),
                 rgb(0.278012, 0.180367, 0.486697),
                 rgb(0.277134, 0.185228, 0.489898),
                 rgb(0.276194, 0.190074, 0.493001),
                 rgb(0.275191, 0.194905, 0.496005),
                 rgb(0.274128, 0.199721, 0.498911),
                 rgb(0.273006, 0.204520, 0.501721),
                 rgb(0.271828, 0.209303, 0.504434),
                 rgb(0.270595, 0.214069, 0.507052),
                 rgb(0.269308, 0.218818, 0.509577),
                 rgb(0.267968, 0.223549, 0.512008),
                 rgb(0.266580, 0.228262, 0.514349),
                 rgb(0.265145, 0.232956, 0.516599),
                 rgb(0.263663, 0.237631, 0.518762),
                 rgb(0.262138, 0.242286, 0.520837),
                 rgb(0.260571, 0.246922, 0.522828),
                 rgb(0.258965, 0.251537, 0.524736),
                 rgb(0.257322, 0.256130, 0.526563),
                 rgb(0.255645, 0.260703, 0.528312),
                 rgb(0.253935, 0.265254, 0.529983),
                 rgb(0.252194, 0.269783, 0.531579),
                 rgb(0.250425, 0.274290, 0.533103),
                 rgb(0.248629, 0.278775, 0.534556),
                 rgb(0.246811, 0.283237, 0.535941),
                 rgb(0.244972, 0.287675, 0.537260),
                 rgb(0.243113, 0.292092, 0.538516),
                 rgb(0.241237, 0.296485, 0.539709),
                 rgb(0.239346, 0.300855, 0.540844),
                 rgb(0.237441, 0.305202, 0.541921),
                 rgb(0.235526, 0.309527, 0.542944),
                 rgb(0.233603, 0.313828, 0.543914),
                 rgb(0.231674, 0.318106, 0.544834),
                 rgb(0.229739, 0.322361, 0.545706),
                 rgb(0.227802, 0.326594, 0.546532),
                 rgb(0.225863, 0.330805, 0.547314),
                 rgb(0.223925, 0.334994, 0.548053),
                 rgb(0.221989, 0.339161, 0.548752),
                 rgb(0.220057, 0.343307, 0.549413),
                 rgb(0.218130, 0.347432, 0.550038),
                 rgb(0.216210, 0.351535, 0.550627),
                 rgb(0.214298, 0.355619, 0.551184),
                 rgb(0.212395, 0.359683, 0.551710),
                 rgb(0.210503, 0.363727, 0.552206),
                 rgb(0.208623, 0.367752, 0.552675),
                 rgb(0.206756, 0.371758, 0.553117),
                 rgb(0.204903, 0.375746, 0.553533),
                 rgb(0.203063, 0.379716, 0.553925),
                 rgb(0.201239, 0.383670, 0.554294),
                 rgb(0.199430, 0.387607, 0.554642),
                 rgb(0.197636, 0.391528, 0.554969),
                 rgb(0.195860, 0.395433, 0.555276),
                 rgb(0.194100, 0.399323, 0.555565),
                 rgb(0.192357, 0.403199, 0.555836),
                 rgb(0.190631, 0.407061, 0.556089),
                 rgb(0.188923, 0.410910, 0.556326),
                 rgb(0.187231, 0.414746, 0.556547),
                 rgb(0.185556, 0.418570, 0.556753),
                 rgb(0.183898, 0.422383, 0.556944),
                 rgb(0.182256, 0.426184, 0.557120),
                 rgb(0.180629, 0.429975, 0.557282),
                 rgb(0.179019, 0.433756, 0.557430),
                 rgb(0.177423, 0.437527, 0.557565),
                 rgb(0.175841, 0.441290, 0.557685),
                 rgb(0.174274, 0.445044, 0.557792),
                 rgb(0.172719, 0.448791, 0.557885),
                 rgb(0.171176, 0.452530, 0.557965),
                 rgb(0.169646, 0.456262, 0.558030),
                 rgb(0.168126, 0.459988, 0.558082),
                 rgb(0.166617, 0.463708, 0.558119),
                 rgb(0.165117, 0.467423, 0.558141),
                 rgb(0.163625, 0.471133, 0.558148),
                 rgb(0.162142, 0.474838, 0.558140),
                 rgb(0.160665, 0.478540, 0.558115),
                 rgb(0.159194, 0.482237, 0.558073),
                 rgb(0.157729, 0.485932, 0.558013),
                 rgb(0.156270, 0.489624, 0.557936),
                 rgb(0.154815, 0.493313, 0.557840),
                 rgb(0.153364, 0.497000, 0.557724),
                 rgb(0.151918, 0.500685, 0.557587),
                 rgb(0.150476, 0.504369, 0.557430),
                 rgb(0.149039, 0.508051, 0.557250),
                 rgb(0.147607, 0.511733, 0.557049),
                 rgb(0.146180, 0.515413, 0.556823),
                 rgb(0.144759, 0.519093, 0.556572),
                 rgb(0.143343, 0.522773, 0.556295),
                 rgb(0.141935, 0.526453, 0.555991),
                 rgb(0.140536, 0.530132, 0.555659),
                 rgb(0.139147, 0.533812, 0.555298),
                 rgb(0.137770, 0.537492, 0.554906),
                 rgb(0.136408, 0.541173, 0.554483),
                 rgb(0.135066, 0.544853, 0.554029),
                 rgb(0.133743, 0.548535, 0.553541),
                 rgb(0.132444, 0.552216, 0.553018),
                 rgb(0.131172, 0.555899, 0.552459),
                 rgb(0.129933, 0.559582, 0.551864),
                 rgb(0.128729, 0.563265, 0.551229),
                 rgb(0.127568, 0.566949, 0.550556),
                 rgb(0.126453, 0.570633, 0.549841),
                 rgb(0.125394, 0.574318, 0.549086),
                 rgb(0.124395, 0.578002, 0.548287),
                 rgb(0.123463, 0.581687, 0.547445),
                 rgb(0.122606, 0.585371, 0.546557),
                 rgb(0.121831, 0.589055, 0.545623),
                 rgb(0.121148, 0.592739, 0.544641),
                 rgb(0.120565, 0.596422, 0.543611),
                 rgb(0.120092, 0.600104, 0.542530),
                 rgb(0.119738, 0.603785, 0.541400),
                 rgb(0.119512, 0.607464, 0.540218),
                 rgb(0.119423, 0.611141, 0.538982),
                 rgb(0.119483, 0.614817, 0.537692),
                 rgb(0.119699, 0.618490, 0.536347),
                 rgb(0.120081, 0.622161, 0.534946),
                 rgb(0.120638, 0.625828, 0.533488),
                 rgb(0.121380, 0.629492, 0.531973),
                 rgb(0.122312, 0.633153, 0.530398),
                 rgb(0.123444, 0.636809, 0.528763),
                 rgb(0.124780, 0.640461, 0.527068),
                 rgb(0.126326, 0.644107, 0.525311),
                 rgb(0.128087, 0.647749, 0.523491),
                 rgb(0.130067, 0.651384, 0.521608),
                 rgb(0.132268, 0.655014, 0.519661),
                 rgb(0.134692, 0.658636, 0.517649),
                 rgb(0.137339, 0.662252, 0.515571),
                 rgb(0.140210, 0.665859, 0.513427),
                 rgb(0.143303, 0.669459, 0.511215),
                 rgb(0.146616, 0.673050, 0.508936),
                 rgb(0.150148, 0.676631, 0.506589),
                 rgb(0.153894, 0.680203, 0.504172),
                 rgb(0.157851, 0.683765, 0.501686),
                 rgb(0.162016, 0.687316, 0.499129),
                 rgb(0.166383, 0.690856, 0.496502),
                 rgb(0.170948, 0.694384, 0.493803),
                 rgb(0.175707, 0.697900, 0.491033),
                 rgb(0.180653, 0.701402, 0.488189),
                 rgb(0.185783, 0.704891, 0.485273),
                 rgb(0.191090, 0.708366, 0.482284),
                 rgb(0.196571, 0.711827, 0.479221),
                 rgb(0.202219, 0.715272, 0.476084),
                 rgb(0.208030, 0.718701, 0.472873),
                 rgb(0.214000, 0.722114, 0.469588),
                 rgb(0.220124, 0.725509, 0.466226),
                 rgb(0.226397, 0.728888, 0.462789),
                 rgb(0.232815, 0.732247, 0.459277),
                 rgb(0.239374, 0.735588, 0.455688),
                 rgb(0.246070, 0.738910, 0.452024),
                 rgb(0.252899, 0.742211, 0.448284),
                 rgb(0.259857, 0.745492, 0.444467),
                 rgb(0.266941, 0.748751, 0.440573),
                 rgb(0.274149, 0.751988, 0.436601),
                 rgb(0.281477, 0.755203, 0.432552),
                 rgb(0.288921, 0.758394, 0.428426),
                 rgb(0.296479, 0.761561, 0.424223),
                 rgb(0.304148, 0.764704, 0.419943),
                 rgb(0.311925, 0.767822, 0.415586),
                 rgb(0.319809, 0.770914, 0.411152),
                 rgb(0.327796, 0.773980, 0.406640),
                 rgb(0.335885, 0.777018, 0.402049),
                 rgb(0.344074, 0.780029, 0.397381),
                 rgb(0.352360, 0.783011, 0.392636),
                 rgb(0.360741, 0.785964, 0.387814),
                 rgb(0.369214, 0.788888, 0.382914),
                 rgb(0.377779, 0.791781, 0.377939),
                 rgb(0.386433, 0.794644, 0.372886),
                 rgb(0.395174, 0.797475, 0.367757),
                 rgb(0.404001, 0.800275, 0.362552),
                 rgb(0.412913, 0.803041, 0.357269),
                 rgb(0.421908, 0.805774, 0.351910),
                 rgb(0.430983, 0.808473, 0.346476),
                 rgb(0.440137, 0.811138, 0.340967),
                 rgb(0.449368, 0.813768, 0.335384),
                 rgb(0.458674, 0.816363, 0.329727),
                 rgb(0.468053, 0.818921, 0.323998),
                 rgb(0.477504, 0.821444, 0.318195),
                 rgb(0.487026, 0.823929, 0.312321),
                 rgb(0.496615, 0.826376, 0.306377),
                 rgb(0.506271, 0.828786, 0.300362),
                 rgb(0.515992, 0.831158, 0.294279),
                 rgb(0.525776, 0.833491, 0.288127),
                 rgb(0.535621, 0.835785, 0.281908),
                 rgb(0.545524, 0.838039, 0.275626),
                 rgb(0.555484, 0.840254, 0.269281),
                 rgb(0.565498, 0.842430, 0.262877),
                 rgb(0.575563, 0.844566, 0.256415),
                 rgb(0.585678, 0.846661, 0.249897),
                 rgb(0.595839, 0.848717, 0.243329),
                 rgb(0.606045, 0.850733, 0.236712),
                 rgb(0.616293, 0.852709, 0.230052),
                 rgb(0.626579, 0.854645, 0.223353),
                 rgb(0.636902, 0.856542, 0.216620),
                 rgb(0.647257, 0.858400, 0.209861),
                 rgb(0.657642, 0.860219, 0.203082),
                 rgb(0.668054, 0.861999, 0.196293),
                 rgb(0.678489, 0.863742, 0.189503),
                 rgb(0.688944, 0.865448, 0.182725),
                 rgb(0.699415, 0.867117, 0.175971),
                 rgb(0.709898, 0.868751, 0.169257),
                 rgb(0.720391, 0.870350, 0.162603),
                 rgb(0.730889, 0.871916, 0.156029),
                 rgb(0.741388, 0.873449, 0.149561),
                 rgb(0.751884, 0.874951, 0.143228),
                 rgb(0.762373, 0.876424, 0.137064),
                 rgb(0.772852, 0.877868, 0.131109),
                 rgb(0.783315, 0.879285, 0.125405),
                 rgb(0.793760, 0.880678, 0.120005),
                 rgb(0.804182, 0.882046, 0.114965),
                 rgb(0.814576, 0.883393, 0.110347),
                 rgb(0.824940, 0.884720, 0.106217),
                 rgb(0.835270, 0.886029, 0.102646),
                 rgb(0.845561, 0.887322, 0.099702),
                 rgb(0.855810, 0.888601, 0.097452),
                 rgb(0.866013, 0.889868, 0.095953),
                 rgb(0.876168, 0.891125, 0.095250),
                 rgb(0.886271, 0.892374, 0.095374),
                 rgb(0.896320, 0.893616, 0.096335),
                 rgb(0.906311, 0.894855, 0.098125),
                 rgb(0.916242, 0.896091, 0.100717),
                 rgb(0.926106, 0.897330, 0.104071),
                 rgb(0.935904, 0.898570, 0.108131),
                 rgb(0.945636, 0.899815, 0.112838),
                 rgb(0.955300, 0.901065, 0.118128),
                 rgb(0.964894, 0.902323, 0.123941),
                 rgb(0.974417, 0.903590, 0.130215),
                 rgb(0.983868, 0.904867, 0.136897),
                 rgb(0.993248, 0.906157, 0.143936))

color.brewer <- function(n)
{
    library(RColorBrewer)
    if(n < 9)
    {
        return (brewer.pal(n, "Set1"))
    }else if(n < 18)
    {
        return (c(brewer.pal(9, "Set1"), brewer.pal(19-n, "Set2")))
    }else if(n < 27)
    {
        return (c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), c(brewer.pal(26-n, "Dark2"))))
    }
}
