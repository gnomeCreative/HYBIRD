

% particl diameter
d=0.002;

% simulation box dimensions
box_width=1.2*d;
box_length=0.035;

% single particle volume (sphere)
v_s=4.0/3.0*pi*(d/2)^2;

% tot particles
p_tot=1000;

% dispersity parameter (needed in order to have particles slightly different in size, and avoid crystallization)
disp=0.05*d;

% largest possible particle diameter (to compute spacing)
d_largest=d+disp;

% spacign a bit larger than max diameters
spacing=d_largest*1.05;

%grid of regularly shaped particles
p_l=floor(box_length/spacing);
p_h=ceil(p_tot/p_l);

% effective number of particles. Should be similar to the one declared
% initially
p_tot=p_h*p_l;

% create particle file and populate it with particles
particleFile='particle_init.dat';
fileID=fopen(particleFile, 'w');

fprintf( fileID,'%i\n',p_tot);

%particle index
index=0;
for i=1:p_l
    for j=1:p_h
        r=0.5*random('unif',d-disp,d+disp);
        x=(i+random('unif',-0.01,0.01)-0.5)*spacing;
        y=box_width/2;
        z=(j+random('unif',-0.01,0.01)-0.5)*spacing;
        fprintf(fileID,'%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',...
            index, 1, r, x, y, z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,...
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        index=index+1;


    end
end

fclose(fileID);






