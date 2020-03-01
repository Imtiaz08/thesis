function genE5=genGoldCode(feedBackReg1,feedBackReg2,loadReg2)
% feedBackReg1 is oct
% feedBackReg2 is oct
% loadReg2 is oct


L=10230;


%% generate



%%

reg1_fb_bin = fliplr(dec2bin(oct2decS(feedBackReg1),15));
reg2_fb_bin = fliplr(dec2bin(oct2decS(feedBackReg2),15));
reg1_fb_vect=reg1_fb_bin-'0';
reg2_fb_vect=reg2_fb_bin-'0';

reg2_start_val_bin = fliplr(dec2bin(oct2decS(loadReg2),14));
reg2_start_val_vect=reg2_start_val_bin-'0';


goldseq = comm.GoldSequence('FirstPolynomial',reg1_fb_vect,...
    'SecondPolynomial',reg2_fb_vect,...
    'FirstInitialConditions', [1 1 1 1 1 1 1 1 1 1 1 1 1 1],...
    'SecondInitialConditions',reg2_start_val_vect,...
    'Index',0,'SamplesPerFrame',L);



genE5 = (step(goldseq))*2-1;



ref = genE5(1:24,:);
% disp('initial bit is');
% disp(ref);


