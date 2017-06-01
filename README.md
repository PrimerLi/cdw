# cdw
A C++ program for solving the mean field equation in Holstein model and periodic Anderson model, with the phonon displacement as the order parameter. 

The Holstein model, which is the simplest model that captures the electron-phonon interaction, is solved using the mean field theory 
method. In the mean field, I have assumed that there is a charge density wave (CDW) order parameter, and 
by imposing the condition that the free energy as a function of the CDW order parameter should be minimized in a real physical system, 
I can obtain a self-consistent equation for the order parameter. With the solution of the self-consistent equation, 
I can find the temperature dependence of the order parameter in different dimensions. What I found is that the larger the dimension, 
the smaller the order parameter, which is consistent with the physical intuition that the CDW is easier to form in lower dimensions, since the perfect nesting condition for the Fermi surface is the most obvious in the one-dimensional system. 

In the periodic  Anderson model with Holstein phonons, I have made the assumption that the Hubbard interaction is negligible, and thus solve the mean field equations for different values of hybridization strengths. What I found is that the electron-phonon interaction is in strong competition with the hybridization effect. This competition is responsible for the many interesting physics that we can see in this model. Because of this competition, the CDW order parameter gets smaller as the hybridization strength gets stronger. 

A detailed description of the model and the methods can be found in the file mft.pdf. 
