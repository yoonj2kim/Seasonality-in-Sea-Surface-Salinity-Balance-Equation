% Load raw data
ncdisp('ArgoSIOForcingAdvectionVerticalMixing.nc')
lat = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','lat');
lon = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','lon');
time = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','time');
sss = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','sio_salinity');
adv = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','sio_advection');
ent = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','sio_vertical_mixing');
forc = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','forcing');
fae = ncread('ArgoSIOForcingAdvectionVerticalMixing.nc','forcing_minus_advection_and_vertical_mixing');

setup_paths
addpath 'C:\Users\yoonj\OneDrive - The Ohio State University\fdaM'

% data missing Jan/Feb 2014 & July 2019
time = time(54:117);
sss = sss(:,:,54:117);
adv = adv(:,:,54:117);
ent = ent(:,:,54:117);
forc = forc(:,:,54:117);
fae = fae(:,:,54:117);

% setting for spline interpolation
m = 200;
t_eval = linspace(0,1,m);
y_num = 4;
t_obs_sss = linspace(0,1,12*y_num+1);

rng      = [0,1];
norder   = 4;
Lfdobj   = int2Lfd(1);
lambda   = 1e-4;

knots_sss    = t_obs_sss;
nbasis_sss   = length(knots_sss) + norder - 2;
hgtbasis_sss = create_bspline_basis(rng, nbasis_sss, norder, knots_sss);
hgtfdPar_sss = fdPar(hgtbasis_sss, Lfdobj, lambda);

%% Register FMAV function to SSS function

latitude = [33 11 -15 48 30 10 45 -25 -25 -20 -8];
longitude = [-150 -150 -170 -158 -30 -30 180 -100 -15 -5 50];

fdobj_sss = cell(length(latitude),5,5);
sss_smooth = cell(length(latitude),5,5);
fdobj_fae = cell(length(latitude),5,5);
fae_smooth = cell(length(latitude),5,5);

rescale_sss = cell(length(latitude),5,5);
rescale_fae = cell(length(latitude),5,5);

gam_com = cell(length(latitude),5,5);
align_dist = zeros(length(latitude),5,5);

sss_scale = zeros(length(latitude),5,5);
fae_scale = zeros(length(latitude),5,5);

for k = 1:length(latitude)
    lat_idx = find(lat == latitude(k));
    lon_idx = find(lon == longitude(k));

        for ii = 1:5
            for jj = 1:5
                sss_temp = squeeze(sss(lon_idx,lat_idx,end-(12*y_num)-2:end-2));
                fae_temp = squeeze(fae(lon_idx,lat_idx,end-(12*y_num)+ii-5:end+jj-5));
                if sum(isempty(sss_temp)) == 0 && sum(isnan(fae_temp)) == 0
                    [fdobj_sss{k},~] = smooth_basis(t_obs_sss, sss_temp, hgtfdPar_sss);
                    sss_smooth{k} = eval_fd(t_eval,fdobj_sss{k});
                    
                    t_obs_fae = linspace(0,1,12*y_num+1+jj-ii);
                    knots_fae    = t_obs_fae;
                    nbasis_fae   = length(knots_fae) + norder - 2;
                    hgtbasis_fae = create_bspline_basis(rng, nbasis_fae, norder, knots_fae);
                    hgtfdPar_fae = fdPar(hgtbasis_fae, Lfdobj, lambda);
                    [fdobj_fae{k,ii,jj},~] = smooth_basis(t_obs_fae, fae_temp, hgtfdPar_fae);
                    fae_smooth{k,ii,jj} = eval_fd(t_eval,fdobj_fae{k,ii,jj});
                    
                    fdot_temp = eval_fd(t_eval,fdobj_sss{k},1);
                    rescale_sss{k} = sss_smooth{k}/trapz(t_eval,abs(fdot_temp));
                    sss_discrete = (sss_temp)/trapz(t_eval,abs(fdot_temp));
                    sss_scale(k,ii,jj) = trapz(t_eval,abs(fdot_temp));

                    fdot_temp = eval_fd(t_eval,fdobj_fae{k,ii,jj},1);
                    rescale_fae{k,ii,jj} = fae_smooth{k,ii,jj}/trapz(t_eval,abs(fdot_temp));
                    fae_discrete = (fae_temp)/trapz(t_eval,abs(fdot_temp));
                    fae_scale(k,ii,jj) = trapz(t_eval,abs(fdot_temp));

                    gam_com{k,ii,jj} = optimum_reparam(f_to_srvf(rescale_fae{k,ii,jj},t_eval),...
                        f_to_srvf(rescale_sss{k},t_eval),t_eval);
                    align_dist(k,ii,jj) = trapz(t_eval,(warp_q_gamma(f_to_srvf(rescale_fae{k,ii,jj},...
                        t_eval),t_eval,gam_com{k,ii,jj})-f_to_srvf(rescale_sss{k},t_eval)).^2);
                else
                    align_dist(k,ii,jj) = NaN;
                end
            end
        end
end

%% Find optimum start and end time point

rescale_fae_min_dist = cell(length(latitude),1);
sss_scale_min_dist = zeros(length(latitude));
fae_scale_min_dist = zeros(length(latitude));

gam_com_min_dist = cell(length(latitude),1);
idx_min_dis = zeros(length(latitude),2);
align_min_dist = zeros(length(latitude));

for k = 1:length(latitude)
        align_dist_temp = squeeze(align_dist(k,:,:));
        align_dist_temp(align_dist_temp==0) = NaN;
        align_min_dist(k) = min(min(align_dist_temp));
        [r,c] = find(squeeze(align_dist(k,:,:))==align_min_dist(k));
        if isempty(r) == 0 && isempty(c) ==0
            gam_com_min_dist{k} = gam_com{k,r,c};
            idx_min_dis(k,:) = [r,c];
            rescale_fae_min_dist{k} = rescale_fae{k,r,c};
            sss_scale_min_dist(k) = sss_scale(k,r,c);
            fae_scale_min_dist(k) = fae_scale(k,r,c);
        end
end

%% figures

t = t_eval;
latitude = [33 11 -15 48 30 10 45 -25 -25 -20 -8];
longitude = [-150 -150 -170 -158 -30 -30 180 -100 -15 -5 50];
phase_diff = [2.9253 3.6001 3.998 2.0997 3.2893 2.6204 1.5624 2.6289 2.7329 3.5056 2.2338];


for kk = 1:length(latitude)
    idx_temp = squeeze(idx_min_dis(kk,:));
    r = idx_temp(1);c = idx_temp(2);
    
    f_1 = sss_scale_min_dist(kk)*rescale_sss{kk};
    f_2 = fae_scale_min_dist(kk)*rescale_fae_min_dist{kk};
    
    t_1 = linspace(0,1,m);
    t_2_l = (r-3)/(y_num*12);
    t_2_u = 1+(c-3)/(y_num*12);
    t_2 = linspace(t_2_l,t_2_u,m);
    
    gam = gam_com_min_dist{kk};
    gam_i = invertGamma(gam);
    K = 20;
    activeIndeces = round(linspace(1,m,K)); 
    cgam = interp1(gam_i,t,t);
    f_1_warp = warp_f_gamma(f_1,t,gam_i);

    ylim_l = [min(f_1)-1.3*(max(f_1)-min(f_1)) max(f_1)+sqrt(var(f_1))];
    ylim_r = [min(f_2)-3*sqrt(var(f_2)) max(f_2)+1.3*(max(f_2)-min(f_2))];

    xlim_l = [min([0,t_2_l]) max([1,t_2_u])];
    
    fig = figure;
    left_color = [0, 0.5, 0];
    right_color = [1.00 0.54 0.00];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    hold off
    yyaxis left
    p1 = plot(t_1,f_1,'LineWidth',2);
    ylim(ylim_l)
    xlim(xlim_l)
    ylabel('SSS','FontSize',12)

    yyaxis right
    p2 = plot(t_2,f_2,'LineWidth',2);
    for k = 1:K
        k_i = activeIndeces(k);
        f_1_proj = ylim_r(1)+(ylim_r(2)-ylim_r(1))*(f_1_warp(k_i)-ylim_l(1))/(ylim_l(2)-ylim_l(1));
        line([cgam(k_i) t_2(k_i)],[f_1_proj f_2(k_i)],'Color',[0.65 0.65 0.65],'LineStyle','--')
    end
    ylim(ylim_r)
    ylabel('FMAV','FontSize',12)
    set(gca,'xtick',[0 1],...
        'xticklabel',{'2015 Apr' '2019 Apr'},'FontSize',12)
    set(gcf,'Units','Inches','Position', [0, 0, 4,2.5], 'PaperUnits','Inches', 'PaperSize', [4,2.5])

    figure
    hold off
    line([0 1],[phase_diff(kk)/(12*y_num) phase_diff(kk)/(12*y_num)],'LineWidth',2,'Color',[0, 0.4470, 0.7410])
    hold on
    plot(t(11:end-10),smooth(gam(11:end-10)-t(11:end-10),3),'Color','k','LineWidth',2)
    max_lag = max(gam(11:end-10)-t(11:end-10));
    min_lag = min(gam(11:end-10)-t(11:end-10));
    line([0 1],[2/(12*y_num) 2/(12*y_num)],'LineStyle','--','Color',[0.65 0.65 0.65])
    line([0 1],[4/(12*y_num) 4/(12*y_num)],'LineStyle','--','Color',[0.65 0.65 0.65])
    line([0 1],[6/(12*y_num) 6/(12*y_num)],'LineStyle','--','Color',[0.65 0.65 0.65])
    line([0 1],[0 0],'LineStyle','--','Color',[0.65 0.65 0.65])
    line([0 1],[-2/(12*y_num) -2/(12*y_num)],'LineStyle','--','Color',[0.65 0.65 0.65])
    set(gca,'xtick',[0 1],...
        'xticklabel',{'2015 Apr' '2019 Apr'},'FontSize',12)
    set(gca,'ytick',[-2/(12*y_num) 0  2/(12*y_num) 4/(12*y_num) 6/(12*y_num)],...
        'yticklabel',["-2","0", "2","4","6"])
    ylabel('Lag (Month)')
    ylim([-3/(12*y_num) 6.5/(12*y_num)])
    set(gcf,'Units','Inches','Position', [0, 0, 4,2.5], 'PaperUnits','Inches', 'PaperSize', [4,2.5])

    fig = figure;
    left_color = [0, 0.5, 0];
    right_color = [1.00 0.54 0.00];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left
    p1 = plot(t,f_1,'LineWidth',2);
    ylim(ylim_l)
    ylabel('SSS','FontSize',12)

    yyaxis right
    p2 = plot(t,warp_f_gamma(f_2,t,gam),'LineWidth',2);
    hold on
    p3 = plot(t+phase_diff(kk)/(12*y_num),f_2-2.5*sqrt(var(f_2)),'LineStyle','-','Color',[0, 0.4470, 0.7410],'LineWidth',1);
    ylim(ylim_r)
    ylabel('FMAV','FontSize',12)
    set(gca,'xtick',[0 1],...
        'xticklabel',{'2015 Apr' '2019 Apr'},'FontSize',12)
    set(gcf,'Units','Inches','Position', [0, 0, 4,2.5], 'PaperUnits','Inches', 'PaperSize', [4,2.5])
end
