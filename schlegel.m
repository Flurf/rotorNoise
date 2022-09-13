function [p_m_rms, db_m] = schlegel(r,a,phi,beta,l,R,sigma,teta,n_harmonica,b,omega,c)
    
    grad2rad = (2*pi) / 360;
    sigma = sigma*grad2rad;
    teta = teta*grad2rad; 
    beta = beta.*grad2rad;
    phi = phi.*grad2rad;
    
    
    
    
    %%%%%%%%%% calcolo p_rms   
    s = zeros(1,length(r),length(phi));
    u = zeros(1,length(r),length(phi));
    
    p_real_harmonic =zeros(1,n_harmonica) ;
    p_imm_harmonic = zeros(1,n_harmonica) ;
    
    p_m_rms = zeros(1,n_harmonica);
    db_m = zeros(1,n_harmonica);
    
    
    for j=1:n_harmonica
        for n=2:length(r)
             for m=2:length(phi)
                s= sqrt( (R^2) + (r(n)^2) - ( 2*r(n)*R*cos(sigma)*cos(teta-phi(m))) ) ;
                u= j*b*( (omega*s/c) + (a*0.5/r(n)) + phi(m) );
                p_real_harmonic(1,j) = p_real_harmonic(1,j) + (l(n,m)/(j*s*s)) * sin((j*a*b*0.5)/r(n)) * ( (j*b*omega*(1/c)*sin(u)) + cos(u)/s ) *   (sin(beta(n,m))*cos(sigma)*sin( phi(m)-teta )  +  cos(beta(n,m))*sin(sigma))  * r(n)  * (phi(m) - phi(m-1))* (r(n) - r(n - 1)); 
                p_imm_harmonic(1,j) = p_imm_harmonic(1,j) + (l(n,m)/(j*s*s)) * sin((j*a*b*0.5)/r(n)) * ( (j*b*omega*(1/c)*cos(u)) - sin(u)/s ) * ( sin(beta(n,m))*cos(sigma)*sin( phi(m)-teta )  +  cos(beta(n,m))*sin(sigma)) * r(n) *  (phi(m) - phi(m-1))* (r(n) - r(n - 1));
             end
        end   
    end
        
    
    for i=1:n_harmonica
        p_m_rms(i)= ((R/ (2*(2^0.5)*pi^2*a) )) *  sqrt( p_imm_harmonic(i)^2 + p_real_harmonic(i)^2 ) ; 
        db_m(i) = 20*log10(p_m_rms(i)/(0.00002));    
    end
end
