using FileIO, Images, DSP, LinearAlgebra, LaTeXStrings

# Cargar imagen de referencias para normalizar
imagen = load("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4723.JPG")
gris = Gray.(imagen)
gris_filtered = imfilter(gris, Kernel.gaussian(5.0))
normalization = maximum(gris_filtered)

## Funciones para leer imagen y hacer el corte
function Itensidad(path, normalization)
    imagen = load(path)
    gris = Gray.(imagen)
    gris_filtered = imfilter(gris, Kernel.gaussian(5.0))
    gris_filtered = gris_filtered ./ normalization
    #gris_normalized = gris_filtered ./ maximum(gris_filtered)
    return Float64.(gris_filtered)
end

function cut(I1, I2, I3, I4, cols, rens)
    I1 = I1[rens,cols] 
    I2 = I2[rens,cols] 
    I3 = I3[rens,cols] 
    I4 = I4[rens,cols] 

    return I1,I2,I3,I4
end    

## Cargar imagenes
I1 =Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4723.JPG", normalization)
I2 = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4724.JPG", normalization)
I3 = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4725.JPG", normalization)
I4 = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4725.JPG", normalization)

I1_obj = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4788.JPG", normalization)
I2_obj = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4789.JPG", normalization)
I3_obj = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4790.JPG", normalization)
I4_obj = Itensidad("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/DSC_4791.JPG", normalization)

## Graficar imagen para definir corte

heatmap(I1_obj)


## Realizar recorte

cols = 2200:3700
rens = 1100:2400

I1,I2,I3,I4 = cut(I1,I2,I3,I4, cols, rens)

I1_obj,I2_obj,I3_obj,I4_obj = cut(I1_obj,I2_obj,I3_obj,I4_obj, cols, rens)

## Verificar recorte
heatmap(ys,xs,I1_obj, aspect_ratio=:equal, xlims=(ys[1], ys[end]), ylims=(xs[1], xs[end]), xlabel = L"x[m]", ylabel = L"y[m]",colorbar_title = L"I/I_0" )
savefig("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Reporte RETO/intensidad_volcan.png")
## Definir distancias con unidades reales

#x_distance = 2.75 
#y_distance = 0.98 #huevo
x_distance = 3.08
y_distance = 1.07
#d = 0.006
d =  0.0127


theta = atan(y_distance/x_distance)

## Contar numero de franjas en el recorte y definir vectores xs y ys

#numero_franjas = 19
#numero_franjas = 19 # estrella cols = 3700:4700, rens = 1200:2300
numero_franjas = 39 # volcan cols = 2200:3700, rens = 1100:2400
# huevocols = 2500:3700 ,rens = 800:1800
#numero_franjas=16

xs = LinRange( 0,(numero_franjas*d*size(I1_obj)[1])/(size(I1_obj)[2] ) ,size(I1_obj)[1])
ys = LinRange( 0,numero_franjas*d ,size(I1_obj)[2])


## Calcular fases envueltas  y aplicar algoritmo de desenvolvimiento

phi_ref = atan.(- (I2 .- I4) ./ (I1 .- I3))
phi_obj = atan.(- (I2_obj .- I4_obj) ./ (I1_obj .- I3_obj))


phi_unwrap_obj =  unwrap(phi_obj, dims=1:2, range=pi)
phi_unwrap_ref =  unwrap(phi_ref, dims=1:2, range=pi)

phi_unwrap = phi_unwrap_obj .- phi_unwrap_ref

## Graficar fases

gr()
heatmap(ys, xs, phi_unwrap, aspect_ratio=:equal, xlims=(ys[1], ys[end]), ylims=(xs[1], xs[end]), xlabel = L"x[m]", ylabel = L"y[m]")
#savefig("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Reporte RETO/unwrap_volcan.pdf")

## Calcular la topografía

z = (phi_unwrap ./ 2pi) .* (d ./ sin(theta)) 
 


## Graficar la topografía

z_min = minimum(z)
z_max = maximum(z)




gr()
surface(ys,xs,z, camera=(35, 65), grid = false, clims=(0,0.4), aspect_ratio=:equal)
xlims!(0, 0.4)
ylims!(0, 0.4)
zlims!(0, 0.4)


##
plotlyjs()
surface(ys,xs,z, camera=(35, 65), grid = false)
plot3d!([0,0,0],[0,0,0],[0,0,0.4], color = "white", alpha = 0,label = false)
plot3d!([0.4,0,0],[0,0,0],[0,0,0], color = "white", alpha=0,label = false)
plot3d!([0,0,0],[0,0.4,0],[0,0,0], color = "white", alpha=0,label = false)
#savefig("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Fotos exp/volcan.pdf")

##
using GLMakie

#surface(ys,xs,z, colormap=:lajolla ) 


##
fig = Figure()
scene = LScene(fig[1,1], show_axis=false)
plt = surface!(scene, ys, xs, z, colormap=:lajolla)
Colorbar(fig[1,2], plt) 
cb = Colorbar(fig[1,2], plt, label=L"z\ [m]")  
display(fig)
# Apply the update


#save("/Users/emilianoreyna/Downloads/ITESM/Laboratorio Óptica/RETO/Reporte RETO/volcan.png", fig)
##

