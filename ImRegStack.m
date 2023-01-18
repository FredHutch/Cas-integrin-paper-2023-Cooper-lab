function [XYTCStackReg, tformall]=ImRegStack(XYTCStack)

% INPUT:    -XYTC_Stack: 3 or 4D image stack
%OUTPUT:    -XYTC_Stack: 3 or 4D registered stack
%           -tformall: cell array containing tforms
tic


sS=size(XYTCStack);
nx=sS(1);
ny=sS(2);
nz=sS(3);
tformall=cell(nz,1);

XYTCStackReg=zeros(sS,'like',XYTCStack);
XYTCStackReg(:,:,1,1)=XYTCStack(:,:,1,1);
sref=[sS(1) sS(2)];

[o,m]=imregconfig('monomodal');
progressText(0,'Image Registration');
for j=2:nz
    ImFixed=XYTCStackReg(:,:,j-1,1);
    ImMoving=XYTCStack(:,:,j,1);
    tform=imregtform(ImMoving,ImFixed,'translation',o,m);
    tformall{j}=tform;
    ImMovingReg=imwarp(ImMoving,tform,'OutputView',imref2d(sref));
    XYTCStackReg(:,:,j,1)=ImMovingReg;
    if numel(sS)==4
        nC=sS(4);
        for k=2:nC
            ImMoving=XYTCStack(:,:,j,k);
            ImMovingReg=imwarp(ImMoving,tform,'OutputView',imref2d(sref));
            XYTCStackReg(:,:,j,k)=ImMovingReg;
        end
    end
    progressText(j/nz);
end

crop=min(XYTCStackReg(:,:,:,1),[],3);
xminp=sum(crop,1);
xtf=xminp==0;
yminp=sum(crop,2);
ytf=yminp==0;
XYTCStackReg(:,xtf,:,:)=[];
XYTCStackReg(ytf,:,:,:)=[];

toc
