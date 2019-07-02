clear all

% Gaussian beam
%   theta = 40;
  rnum  = 1000; 
  OMEGA = rnum*2*pi;
  rk = OMEGA;
  zi = sqrt(-1);  
% 
% beam waist radius
  w0 = 0.02;  
  [X,Y,Z] = meshgrid(0:0.01:1,0:0.01:1, 0:0.01:1);
  
  
  % Shift the origin
  xprim = X-0.01;
  yprim = Y-0.01;
  zprim = Z-0.01;
  
% Rotation
theta(1) = 35; theta(2)= -45; theta(3) =0;
rad   = theta*pi/180; 
% rad(1) = pi/2-atan(sqrt(3)/3);
% rad(2) = asin(sqrt(3)/3);
% rad(3) = acos(sqrt(3)/3);
sina(1:3)=sin(rad(1:3)); cosa(1:3) = cos(rad(1:3));

Rx = [1 0 0 ; 0 cosa(1) -sina(1) ; 0 sina(1) cosa(1)];
Ry = [cosa(2) 0 sina(2) ; 0 1 0 ; -sina(2) 0 cosa(2)];
Rz = [cosa(3) -sina(3) 0 ; sina(3) cosa(3) 0 ; 0 0 1];

R = Rx*Ry*Rz;


   
    x = R(1,1)*xprim+R(1,2)*yprim+R(1,3)*zprim;
    y = R(2,1)*xprim+R(2,2)*yprim+R(2,3)*zprim;
    z = R(3,1)*xprim+R(3,2)*yprim+R(3,3)*zprim;

%   x = xprim;
%   y = yprim;
%   z = zprim;
  r = sqrt(x.^2+y.^2);
  
  % wavelength
  rl = 2*pi/rk;

% 
% function w
  fact = rl/pi/w0.^2;
  aux = 1 + (fact.*z).^2;
  w = w0*aux.^(0.5);
%     
% phase function
  phi0 = atan(fact.*z);
%    
% beam radius
  zr = z + 1./z./fact.^2;
% 
% pressure
  ze = -r.^2./w.^2 - zi*rk.*z - zi*pi.*r.^2./rl./zr + zi*phi0./2;
  pf = (2/pi./w.^2).^0.25;
  zp = pf.*exp(ze);
  
  zp = real(zp);
  
  for i = 1:size(zp,1)
      for j = 1:size(zp,2)
          for k = 1:size(zp,3)
              if(abs(zp(i,j,k)) < 0.1)
                  zp(i,j,k) = NaN;
              end
          end
      end
  end
  
  
figure
scatter3( X(:), Y(:), Z(:), 3, zp(:))
colorbar
colormap(hsv)
caxis([-2,2])
axis([0 1 0 1 0 1 ])
xlabel('x');ylabel('y');zlabel('z');
az = -12;
el = 20;
view(az,el)





