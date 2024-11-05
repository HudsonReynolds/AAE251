mDiff = 0;

mI1 = 0;

mI2 = 0;

payload = 22800;

mPay = 100;

while mI1 < payload * 0.9

mPay = mPay + mDiff / 10;

fInert1 = 0.12;

fInert2 = 0.08;

dV1 = 3506;

dV2 = 3279;

C1 = 430 * 9.81;

C2 = 315 * 9.81;

MR1 = exp(dV1/C1);

MR2 = exp(dV2/C2);

% stage 2:

mProp2 = mPay * (MR2-1) * (1 - fInert2) / (1 - fInert2 * MR2);

mF = mProp2 / (MR2-1);
 
mI2 = mF + mProp2;
 
m_inert = mF - mPay;

fInert2 = m_inert / (m_inert + mProp2)


% stage 1


mI1 = mPay * ((MR1*(1-fInert1))/(1-fInert1*MR1)) * ((MR2*(1-fInert2))/(1-fInert2*MR2));

mF1 = mI1 / MR1;

mProp1 = mI1 - mF1;

mInert1 = mF1 - mI2;

fInert1 = mInert1 / (mInert1 + mProp1);

mDiff = payload * .9 - mI1;

end

fprintf("Final Payload Mass is %.2f and Initial Mass is %.2f. Payload limit is %.2f\n", mPay, mI1, payload * .9)








