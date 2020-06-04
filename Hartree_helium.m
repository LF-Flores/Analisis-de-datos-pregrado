% Base elegida para comenzar el algoritmo
alpha = [0.297104; 1.236745; 5.749982; 38.216677];
N = length(alpha);
tic()
% Construye S T A H
for p=1:N
  for q=1:N
    S(p,q) = (pi/(alpha(p)+alpha(q)))^(3/2);
    T(p,q) = 3*alpha(p)*alpha(q)*pi^(1.5)/(alpha(p)+alpha(q))^(5/2);
    A(p,q) = -4*pi/(alpha(p)+alpha(q));
    H(p,q) = T(p,q)+A(p,q);
  end
end

% Construye Q
for p=1:N
  for q=1:N
    for r=1:N
      for s=1:N
        Q(p,q,r,s) = 2*pi^(2.5)/(((alpha(p)+alpha(q))*(alpha(r)+alpha(s))*(alpha(p)+alpha(q)+alpha(r)+alpha(s))^(0.5)));
      end
    end
   end  
end

%% autoconsistencia de Hartree-Fock 
c = zeros(4,1); c(1) = 1; % Vector de integración inicializado
c1 = zeros(4,1); % Actualizador de C
F = zeros(4,4); v = zeros(4,4);
testing = 0;
const = 1; NormS= 0;

while(const >=1e-10)
  NormS = c'*S*c;
  % Normalización
  c = c/((NormS)^(0.5));

  % Crea la matriz F de Fock
    for p=1:N
      for q=1:N
        F(p,q) = H(p,q);
        for s=1:N
          for r=1:N
            F(p,q) = F(p,q) + c(r)*c(s)*Q(p,q,r,s);
          end
        end
      end
    end

  % Encontrando un autovector para el min(espectro)
  c1 = c;    % Actualizando el autovector
  % F*V = S*V*D 
  % V = cv, D = E
  [V, D] = eig(F,S);
  c = V(:,1);

  % Test de convergencia
  const = (c-c1)' *(c-c1);
end  


E = 0;
% Encontrando la energía del estado base del átomo de hidrógeno (en unidades atómicas)
% Esta debe ser igual a -2.85516 a.u (verificado)
for p=1:N
  for q=1:N
    E = E + 2*c(p)*c(q)*H(p,q);
    for r=1:N
      for s=1:N
        E = E + Q(p,q,r,s)*c(p)*c(q)*c(r)*c(s);
      end
    end
  end
end
toc()
