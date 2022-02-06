function [t0, x0, u0] = shiftLYC(t0, x0, u, R)

st = x0;
con = u(:,1);

st_next = (st(1:R.Nl-1)+R.T*st(R.Nl+3:2*R.Nl+1));
st_next = [st_next;st(R.Nl)+R.T*st(2*R.Nl+2)];
st_next = [st_next;st(R.Nl+1)+R.T*(st(2*R.Nl+3)*cos(st(R.Nl))-st(2*R.Nl+4)*sin(st(R.Nl)))];
st_next = [st_next;st(R.Nl+2)+R.T*(st(2*R.Nl+3)*sin(st(R.Nl))+st(2*R.Nl+4)*cos(st(R.Nl)))];
%     st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*(R.cp/R.m*st(2*R.Nl+3)*R.A*R.D'*st(1:R.Nl-1)+1/R.m*R.D*R.D'*con)];
st_next = [st_next;st(R.Nl+3:2*R.Nl+1)+R.T*con];
st_next = [st_next;st(2*R.Nl+2)+R.T*(-R.lambda1*st(2*R.Nl+2)+R.lambda2/(R.Nl-1)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];
st_next = [st_next;st(2*R.Nl+3)+R.T*(-R.ct/R.m*st(2*R.Nl+3)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+4)*R.e_e'*st(1:R.Nl-1)...
    -R.cp/(R.Nl*R.m)*st(1:R.Nl-1)'*R.A*R.D_D*st(R.Nl+3:2*R.Nl+1))];
st_next = [st_next;st(2*R.Nl+4)+R.T*(-R.cn/R.m*st(2*R.Nl+4)+2*R.cp/(R.Nl*R.m)*st(2*R.Nl+3)*R.e_e'*st(1:R.Nl-1))];

x0 = full(st_next);
t0 = t0 + R.T;
%׼����һ�����Ƶ����ſ��ƣ���Ϊu[:, 0]�Ѿ����ɣ��ͼ򵥵ذѺ���Ľ����ǰ
u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end