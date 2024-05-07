function [Stress, Stiffness]=GetHDStressStiffness(E0,SigmaMax,Strain)
%Caculate stress and stiffness from the Modified Hardin Drenvich Stress Strain relationship   
Stress= E0*Strain/(1+E0*Strain/SigmaMax);
Stiffness=E0*SigmaMax^2/(E0*Strain+SigmaMax)^2;  
end