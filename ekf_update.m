function [X, Cx, a, Ca] = ekf_update(X, Z, Cx, D, a, Ca, Q, F, H)

% Prediction for state vector and covariance:
X = F*X;
Cx = F * Cx * F.' + Q;

avec = [a, (D(1,1))^(1/2), (D(2,2))^(1/2)].';
Ca = eye(3) * Ca * eye(3)';

if sum(sum(~isnan(Z))) > 0  
    for j = 1:sum(~isnan(Z(1,:)),2)

        %R = calcrotmat(X);
        R = [cos(avec(1)), -sin(avec(1)); sin(avec(1)), cos(avec(1))];
        Drot = R * D * R.';
        

        Ex = H*X;

        % Compute conventional Kalman covariance matrices:
        Cxz = Cx*H.';
        Czz = H*Cx*H.' + Drot;

        % Correction based on observation:
        X = X + Cxz*(Czz)^(-1)*(Z(:,j)-Ex);
        Cx = Cx - Cxz*(Czz)^(-1)*Cxz.';
        Cx = (Cx + Cx.')/2;

        % Compute pseudo-measurement and expectation:
        %Ztilde = kron((Z - X(1:2)).',(Z - X(1:2)).');
        Zshift = Z(:,j) - Ex;
        Ztilde = [eye(2), zeros(2,2); 0 0 0 1]*kron(Zshift,Zshift);
        %Ztilde = [(Z(1,:) - X(1)).^2; (Z(2,:) - X(2)).^2; (Z(1,:) - X(1)) .* (Z(2,:) - X(2))];
        %expect_Ztilde = [Czz(1,1), Czz(2,2), Czz(1,2)].';
        expect_Ztilde = [Czz(1,1); Czz(1,2); Czz(2,2)];

        % Compute extended covariance matrices:
        elem1 = 3*(Czz(1,1))^2; elem2 = 3*(Czz(2,2))^2;
        elem3 = Czz(1,1) * Czz(2,2) + 2*(Czz(1,2))^2; elem4 = 3*Czz(1,1)*Czz(1,2);
        elem5 = 3*Czz(2,2)*Czz(1,2);
        Czztilde = [elem1, elem4, elem3; elem4, elem3, elem5; elem3, elem5, elem2];

        expect_Ja = [-sin(2*avec(1)), (cos(avec(1)))^2, (sin(avec(1)))^2; cos(2*avec(1)), sin(2*avec(1)), -sin(2*avec(1)); sin(2*avec(1)), (sin(avec(1)))^2, (cos(avec(1)))^2];
        expect_Ja = expect_Ja * diag([D(1,1)/100 - D(2,2)/100, 2*(D(1,1))^(1/2)/100, 2*(D(2,2))^(1/2)/100]);
        Caztilde = Ca * expect_Ja.';
        % Correction based on observation:
        avec = avec  + Caztilde*(Czztilde)^(-1)*(Ztilde-expect_Ztilde);
        avec(1) = mod(avec(1), 2*pi);
        Ca = Ca - Caztilde * (Czztilde)^(-1) * Caztilde.';
        Ca = (Ca + Ca')/2;
    end
    a = avec(1);
    %a = mod(avec(1), 2*pi);
end


end