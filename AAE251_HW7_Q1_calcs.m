%% option 2

m_pay = 5200;

f_inert = 0.16;

dV = 4350;

C = 305 * 9.81;

m_prop = m_pay * (exp(dV/C)-1) * (1 - f_inert) / (1 - f_inert * exp(dV/C));

MR = exp(dV/C);

m_f = m_prop / (MR-1);

m_i = m_f + m_prop;

m_inert = m_f - m_pay;

f_inert = m_inert / (m_inert + m_prop);

%% Option 1

%stage 3

m_pay = 5200;

f_inert = 0.12;

dV1 = dV * 0.75;

dV2 = dV * 0.25;

C = 282 * 9.81;

MR3 = exp(dV2/C);

m_prop = m_pay * (MR3-1) * (1 - f_inert) / (1 - f_inert * MR3);

m_f = m_prop / (MR3-1);
 
m_i2 = m_f + m_prop;
 
m_inert = m_f - m_pay;

f_inert = m_inert / (m_inert + m_prop)

% stage 2:

MR2 = exp(dV1/C);

m_i = m_pay * (MR3 * (1-f_inert) * MR2 * (1 - f_inert)) / ((1 - f_inert * MR2) * (1 - f_inert * MR3))

m_f = m_i / MR2;

m_prop = m_i - m_f;

m_inert = m_f - m_i2;

f_inert = m_inert / (m_inert + m_prop)



