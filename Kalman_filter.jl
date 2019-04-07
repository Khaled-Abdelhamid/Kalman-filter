#=
x̲ the state vector x̲⁺ posteriori value of the state vector x̲⁻ priori value of the state vector

xₖ = Axₖ₋₁ + Buₖ + wₖ

zₖ = Cxₖ + vₖ

where :

1- u is the input of the system

2- A is the transition matrix

3- B is the the input matrix

4- w is a random variable folowing a guassian ditribution with covariance Q and 0 mean , which mainly represents the proccess noise

5- z is the observation vector which is linearly dependent of the state vector x

6- C is the observation matrix

7- v is a random variable folowing a guassian ditribution with covariance R and 0 mean ,which mainly represents the observation noise

8- P is the covariance of the state vector estimate

9- K is the kalman gain factor
/////////////////////////////////////////////////////////////////////////////////
the system i want to model is car system in which i get its velocity and position
assuming that it moves with constant velocity

the distance equation is Xₖ=X₀+V*dt
where
1-Xₖ:the current position
2-X₀:the intial position
3-V:the current velocity
4-dt: the time duration
////////////////////////////////////////////////////////////////////////////////
=#
using Plots
gr();


mutable struct kal_Data
          A # the transition matrix
          B # the the input matrix
          C # the observation matrix
          Q # process noise's covariance matrix
          R # measure noise's covariance matrix

          u # the input matrix
          z # the observation vector
          x # the state vector
          P # state vector covariance matrix

             kal_Data() = new()
             end

kal=kal_Data() #creating a new instance

function kalman!(kal_struct::kal_Data)

# This is the code which implements the discrete Kalman filter:
   #Prediction for state vector and covariance:
   kal_struct.x = kal_struct.A * kal_struct.x + kal_struct.B * kal_struct.u;
   kal_struct.P = kal_struct.A * kal_struct.P * kal_struct.A' + kal_struct.Q;


   K =(kal_struct.P * kal_struct.C') * inv( kal_struct.C * kal_struct.P * kal_struct.C' + kal_struct.R); # calculation of the kalman gain factor

   #Correction based on observation:
   kal_struct.x = kal_struct.x + K * (kal_struct.z - kal_struct.C * kal_struct.x);
   kal_struct.P = kal_struct.P - K * kal_struct.C * kal_struct.P;
end

#simulation //////////////////////////////////////////////////
dt=0.05     #the sampling time
duration=15 #simulation duration in seconds
t=0:dt:duration #the time range of the simulation

Vᵣ=5 # this is the velocity used in order to calculate the real position
x₀=0 #the intitial position used

Xᵣ=hcat(x₀.+ Vᵣ.*Array(t),Vᵣ.*ones(Float64,length(t)))# the real value that will be used to simulate the measurment
Xₖ=Array{Float64,2}(undef, length(t),2) #preallocation of the sdtate variables array
kal.x=[0.01,0.01]
Zₖ=Array{Float64,2}(undef, length(t),2)   #preallocation of the observation array
u=zeros(Float64, length(t),2)#preallocation of the input array,we are not giving any forcing input to the system so we make it qual to zero

kal.A= [1 dt ;   #intializing the system dynamics array
        0  1 ]

kal.B= [0  0 ;   #intializing the the input array, here we don't have any forcing intervention over the states so it equals zero
        0  0 ]

kal.P = [1  0 ;   #intializing the uncertitnty martix of the model
         0  1 ]


kal.Q = [0  0 ;   #the covariance matrix for the process noise
         0  0 ]

kal.C=  [1 0;    #the observation matrix
         0 1]

kal.R=[1 0;     # the covariance of the measurment noise
       0 1]


for i in 1:length(t)

   kal.u=u[i,:]
   kal.z = kal.C * Xᵣ[i,:] + kal.R*randn(Float64,length(Xᵣ[1,:]),1)
   Zₖ[i,:]=kal.z
   kalman!(kal)
   Xₖ[i,:]=kal.x

end

#=
plot(Xₖ[:,1],title="Kalman filter",label=["pridection"],lw=.5)
plot!(Zₖ[:,1],label=["measurment"],lw=.5)
plot!(Xᵣ[:,1],label=["Real"],lw=.5)
gui()
=#
#now we will display the measured, real and predicted distances
data = [Xₖ[:,1] Zₖ[:,1] Xᵣ[:,1]]  # We put labels in a row vector: applies to each series
labels = ["pridection" "measurment" "Real" ]
plot(t, data, title="Kalman filter (distamce)" ,label = labels,markersize = 10)
xlabel!("time (s)")
ylabel!("Distance (m)")
savefig( "d:\\figure1")

#now we will display the measured, real and predicted Velocities
data = [Xₖ[:,2] Zₖ[:,2] Xᵣ[:,2]]
labels = ["pridection" "measurment" "Real" ]
plot(t, data, title="Kalman filter velocity" ,label = labels,markersize = .5)
xlabel!("time (s)")
ylabel!("velocity (m/s)")
savefig( "d:\\figure2")
