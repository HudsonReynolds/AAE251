%mass for each launch configuration:
payload = [22800, 63800, 16800, 27200, 21000, 95000];

launchVehicle = ["Falcon 9 LEO", "Falcon Heavy", "Falcon Heavy Venus Transfer", "Vulcan Centaur", "Ariane 5", "SLS"];


for i = 1:length(payload)

    capacity = payload(i) * 0.98;
    
    mPay = 100; %initial payload mass
    
    residuals = 1.015; % tank residuals of 1.5%

    % initialize the mass and mass difference to zero to run the loop:
    mI1 = 0;
    mDiff = 0;

    % define constants
    fInert1 = 0.24;
    fInert2 = 0.1448;
    dV1 = 3506 * residuals;
    dV2 = 3279 * residuals;    
    C1 = 453 * 9.81;    
    C2 = 316 * 9.81;    
    MR1 = exp(dV1/C1);    
    MR2 = exp(dV2/C2);

    % run loop for Falcon Heavy, no first stage, direct transfer to Venus:
    if i == 3

        mI2 = 0;
        mDiff = 0;

        while mI2 < capacity
    
            mPay = mPay + mDiff / 20 + .001;

            % stage 2:
            
            mProp2 = mPay * (MR2-1) * (1 - fInert2) / (1 - fInert2 * MR2);
            
            mF = mProp2 / (MR2-1);
             
            mI2 = mF + mProp2;
             
            m_inert = mF - mPay;
            
            fInert2 = m_inert / (m_inert + mProp2);
            
            mDiff = capacity - mI2;
    
        end
        fprintf("Final Payload Mass for %s is %.2f and Initial Mass is %.2f. Payload limit is %.2f\n",launchVehicle(i), mPay, mI2, capacity);

    % run loop for conditions to low earth orbit:
    else

        mI1 = 0;
        mDiff = 0;

        while mI1 < capacity

            mPay = mPay + mDiff / 20;

            % stage 2:
            
            mProp2 = mPay * (MR2-1) * (1 - fInert2) / (1 - fInert2 * MR2);
            
            mF = mProp2 / (MR2-1);
             
            mI2 = mF + mProp2;
             
            m_inert = mF - mPay;
            
            fInert2 = m_inert / (m_inert + mProp2);
            
            
            % stage 1
            
            
            mI1 = mPay * ((MR1*(1-fInert1))/(1-fInert1*MR1)) * ((MR2*(1-fInert2))/(1-fInert2*MR2));
            
            mF1 = mI1 / MR1;
            
            mProp1 = mI1 - mF1;
            
            mInert1 = mF1 - mI2;
            
            fInert1 = mInert1 / (mInert1 + mProp1);
            
            mDiff = capacity - mI1;
            
        end
    end

    fprintf("Final Payload Mass for %s is %.2f and Initial Mass is %.2f. Payload limit is %.2f\n",launchVehicle(i), mPay, mI1, capacity)

end