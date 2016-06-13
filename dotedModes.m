function [ MS MA ] = dotedModes()
% Результати:
% MS - масив симетричних мод
% MA - масив антисиметричних мод
% Структура масивів МА чи MS:
% MS(Nmodes).wk(1,:) - масив частот w
% MS(Nmodes).wk(2,:) - масив хвильових чисел k
% MS(Nmodes).dw(1,:) - масив початков. значень w без екстраполяції
% MS(Nmodes).dw(2,:) - масив початков. значень w обчислених за екстраполяцією 
% наприклад, MS(2).wk(1,:) - 1-ший рядок заповнений частотами моди А1, оскільки є мода А0 

% Шукає корені характеристичного рівняння для пошуку мод пластинчатих хвиль

global d cl ct;
d = 0.0093;   % товщина пластини, м
%cl = 5960; ct = 3260;   %volume waves velocities
E = 210E9;   nu = 0.30;   ro = 7950;    % 2
[cl ct] = soundVelocity(E, nu, ro);

analysis = menu('Вибір розрахунку', 'Дисперсія пластинчатих хвиль', ...
    'Розподіл  поля пластинчатих хвиль по товщині', ...
    'Дисперсія кільцевих хвиль');
if analysis == 1    %дисперсія в пластині
    % w_min w_max - мінімальна максимальна частота діапазону пошуку
    % k_min k_max - мінімальне максимальне хвильове число діапазону пошуку
    w_min = 5; w_max = 7E6; k_min = 1; k_max = 3000;
    [MA MS] = plateModes(w_min, w_max, k_min, k_max);
    
elseif analysis == 2    % профіль поля по товщині в пластині
    w = 2.005099110906770e+06; 
    ka = 708; ks = 708;
    [MA MS] = plateWaveProfile(w, ks, ka, E, nu);
   
elseif analysis == 3    % SH дисперсія в трубі
    w_min = 5; w_max = 7E6; k_min = 1; k_max = 3000;
    R1 = 0.273;
    [MA MS] = dispertionCircumferencial(R1, w_min, w_max, k_min, k_max);
      
end

 MA.ct = ct; MA.cl = cl;
 MS.ct = ct; MS.cl = cl;

end

%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------

function [MA MS] = dispertionCircumferencial(R1, w_min, w_max, k_min, k_max)
% Головна функція обчислення дисперсії в трубі. МА- SH хвиля;   MS - Lamb Circumferen. хвиля
precision = 1E-6;  % поріг, якщо результат функції є меншим то записуємо
dk = 4;    % крок по осі к
[MA] = iterKPipe('dispersionSH', w_min, w_max, k_min, k_max, dk, precision, R1);
[MS] = iterKPipe('dispersionLamb', w_min, w_max, k_min, k_max, dk, precision, R1);

end

%----------------------------------------------------------------------------------------------

function [MM] = iterKPipe(fun, w_min, w_max, k_min, k_max, dk, precision, R1)
% здійснює пошук коренів по осі к для діапазону частот
% Результати
% MM.w - матриця у кожному рядку всі корені для одного значення к
% MM.k - стовпчик значень к

MM.k=0;
MM.w=0;
index = 0;
MM.k = k_min;
for kk = k_min : dk : k_max
    [ VW ] = rapidSearchPipe(fun, w_min, w_max, kk, precision, R1);
    if ~isempty(VW)
        index = index + 1;
        MM.k(index) = kk;
        for ww=1:length(VW)
            MM.w(index, ww) = VW(ww);
        end
    end  
end
MM.k = MM.k';

end

%----------------------------------------------------------------------------------------------

function [ VW ] = rapidSearchPipe(fun, w_min, w_max, k, precision, R1)
% Шукає перші точки мод на відрізку w_min .. w_max
% fun - функція для пошуку кореня, fun(R1, kth, w).
% w_min, w_max - мінімальна (максимальна) частота діапазону пошуку
% precision - точність пошуку кореня алгоритмом
% R1 - внутрішній радіс труби
% Результати:
% VW - вектор знайдених коренів мод

% пошук по діапазону w_min ... w_max, розбитому на 1000 відрізікв
rootsFound = 0; % кількість знайдених коренів
iterw = 0;
stepw = w_max/1000;
for ww=w_min:stepw:(w_max-stepw)
    f1 = feval(fun, R1, k, ww);
    f2 = feval(fun, R1, k, ww + stepw);
    if (f1 * f2) < 0
            rootsFound = rootsFound+1;
            iterw(1,rootsFound) = ww;
            iterw(2,rootsFound) = ww + stepw;
    end
end

% уточнення w та виокремлення мод
for iterat = 1:5    % 5 разів уточнююється в 10 раз зменшуючи відрізок, тому 10*10*10*10*10=100 000
    for modeii = 1:rootsFound
        modeLeft = iterw(1, modeii);
        modeRight = iterw(2, modeii);
        stepw = (modeRight - modeLeft) / 10;

        for ww = modeLeft:stepw:(modeRight - stepw)
            f1 = feval(fun, R1, k, ww);
            f2 = feval(fun, R1, k, ww + stepw);
            if (f1 * f2) < 0
                iterw(1,modeii) = ww;
                iterw(2,modeii) = ww + stepw;
                break;
            end
        end

    end
end

VW = zeros(1, rootsFound);
% здійснюємо пошук кореня
for ii=1:rootsFound
    funsearch = @(ww) feval(fun, R1, k, ww);
    [ww] = bisect(funsearch, iterw(1,ii), iterw(2,ii), precision);
    VW(ii) = ww;      % перший рядок w

end

end

%----------------------------------------------------------------------------------------------

function res = dispersionSH(R1, kth, w)
% Обчислює дисперсію кільцевих горизонтально-поляризованих хвиль
% R1 - внутрішній радіс труби
% kth - кутове хвильове число
% w - частота
global d ct;
h = 1E-4;   % константа для дифференціювання
kt = w/ct;
R2 = R1 + d;
w3 = kt*R1;
w4 = kt*R2;
Jr1_prime = (besselj(kth , w3 + h) - besselj(kth , w3 - h)) / (2*h);
Yr1_prime = (bessely(kth , w3 + h) - bessely(kth , w3 - h)) / (2*h);
Jr2_prime = (besselj(kth , w4 + h) - besselj(kth , w4 - h)) / (2*h);
Yr2_prime = (bessely(kth , w4 + h) - bessely(kth , w4 - h)) / (2*h);
res = Jr1_prime * Yr2_prime - Jr2_prime * Yr1_prime;
end

%----------------------------------------------------------------------------------------------

function res = dispersionLamb(R1, kth, w)
% Обчислює дисперсію кільцевих горизонтально-поляризованих хвиль
% R1 - внутрішній радіс труби
% kth - кутове хвильове число
% w - частота
global d cl ct;
h = 1E-4;   % константа для дифференціювання
kt = w/ct;
kl = w/cl;
R2 = R1 + d;
w1 = kl*R1;
w2 = kl*R2;
w3 = kt*R1;
w4 = kt*R2;
Jp = @(ww) (besselj(kth , ww + h) - besselj(kth , ww - h)) / (2*h);
Yp = @(ww) (bessely(kth , ww + h) - bessely(kth , ww - h)) / (2*h);
JJ = @(ww) besselj(kth, ww);
YY = @(ww) bessely(kth, ww);
f5 = (2*kth*kth) - (w4*w4);
f6 = (2*kth*kth) - (w3*w3);

ch1 = f5*f5*f6*f6*( (JJ(w2)*YY(w1)) - (JJ(w1)*YY(w2)) ) * ( (JJ(w4)*YY(w3)) - (JJ(w3)*YY(w4)) );
ch2 = 16.0*(kth^4)*w1*w2*w3*w4*( (Jp(w2)*Yp(w1)) - (Jp(w1)*Yp(w2)) ) * ( (Jp(w4)*Yp(w3)) - (Jp(w3)*Yp(w4)) ); 
ch3 = -4*kth*kth*f5*f5*w1*w3*( (JJ(w2)*Yp(w1)) - (Jp(w1)*YY(w2)) ) * ( (JJ(w4)*Yp(w3)) - (Jp(w3)*YY(w4)) );
ch4 = -4*kth*kth*f6*f6*w2*w4*( (Jp(w2)*YY(w1)) - (JJ(w1)*Yp(w2)) ) * ( (Jp(w4)*YY(w3)) - (JJ(w3)*Yp(w4)) );
ch5 = (32.0 / (pi*pi)) * (kth*kth*f5*f6);

res = ch1 + ch2 + ch3 + ch4 + ch5;
end

%----------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------

function [MA MS] = plateWaveProfile(w, ks, ka, E, nu)
% Шукає профіль пластинчатих хвиль по товщині пластини u та Sigma
% plot(imag(MS.sig33)+real(MS.sig33), MS.x);
global d cl ct;

% Константи Ляме
mu = E / (2 * ( 1 + nu ));
lambda = (nu*E) / ((1+nu)*( 1-(2*nu) ));

% Розрахунок змінних p та q
q2=(w/ct)^2-ks^2;   % деякі змінні в квадраті
p2=(w/cl)^2-ks^2;   % деякі змінні в квадраті
h = d/2;    % 2h=d товщина стінки пластини
if ks > (w/ct)
    qq = sqrt(q2)/i; pp = sqrt(p2)/i;
    qq2 = -q2;
elseif ((w/ct)>ks) && (ks>(w/cl))
    qq = sqrt(q2); pp = sqrt(p2)/i;
    qq2 = q2;
else
    qq = sqrt(q2); pp = sqrt(p2);
    qq2 = q2;
end

% Граничні умови для визначення амплітуд
% Симетрична мода
cf1 = -2*mu*ks*pp*i*sin(pp*h);   % для Sigma31
cf2 = mu*(ks^2 - qq2)*sin(qq*h);
A2 = 1;
B1 = - A2 * (cf1/cf2);

% Симетрична мода
cf11 = @(xx) i*ks*A2*cos(pp*xx);  % u1
cf12 = @(xx) qq*B1*cos(qq*xx);
cf13 = @(xx) -pp*A2*sin(pp*xx);  % u3
cf14 = @(xx) -i*ks*B1*sin(qq*xx);
cf15 = @(xx) -mu*2*i*ks*pp*A2*sin(pp*xx); % Sigma31
cf16 = @(xx) mu*(ks^2 - qq2)*B1*sin(qq*xx);
cf17 = @(xx) -lambda*(ks^2 + pp^2)*A2*cos(pp*xx);     % Sigma33
cf18 = @(xx) -2*mu*(pp^2)*A2*cos(pp*xx);
cf19 = @(xx) -2*mu*i*ks*qq*B1*cos(qq*xx);


% Розрахунок змінних p та q
q2a=(w/ct)^2-ka^2;   % деякі змінні в квадраті
p2a=(w/cl)^2-ka^2;   % деякі змінні в квадраті
h = d/2;    % 2h=d товщина стінки пластини
if ka > (w/ct)
    qqa = sqrt(q2a)/i; ppa = sqrt(p2a)/i;
    qq2a = -q2a;
elseif ((w/ct)>ka) && (ka>(w/cl))
    qqa = sqrt(q2a); ppa = sqrt(p2a)/i;
    qq2a = q2a;
else
    qqa = sqrt(q2a); ppa = sqrt(p2a);
    qq2a = q2a;
end

% Антсиметрична мода
cf1 = -2*mu*ka*ppa*i*cos(ppa*h);   % для Sigma31
cf2 = mu*(ka^2 - qq2a)*cos(qqa*h);
A1 = 1;
B2 = - A1 * (cf1/cf2);

% Антисиметрична мода
cf21 = @(xx) i*ka*A1*sin(ppa*xx);  % u1
cf22 = @(xx) -qqa*B2*sin(qqa*xx);
cf23 = @(xx) ppa*A1*cos(ppa*xx);  % u3
cf24 = @(xx) -i*ka*B2*cos(qqa*xx);
cf25 = @(xx) mu*2*i*ka*ppa*A1*cos(ppa*xx); % Sigma31
cf26 = @(xx) mu*(ka^2 - qq2a)*B2*cos(qqa*xx);
cf27 = @(xx) -lambda*(ka^2 + ppa^2)*A1*sin(ppa*xx);     % Sigma33
cf28 = @(xx) -2*mu*(ppa^2)*A1*sin(ppa*xx);
cf29 = @(xx) 2*mu*i*ka*qqa*B2*sin(qqa*xx);

xx1 = 0:(h/100):h;
xx2 = -h:(h/100):0;
xx = cat(2,xx2,xx1);

for ii=1:length(xx)
    MS.u1(ii) = cf11(xx(ii)) + cf12(xx(ii));
    MS.u3(ii) = cf13(xx(ii)) + cf14(xx(ii));
    MS.sig31(ii) = cf15(xx(ii)) + cf16(xx(ii));
    MS.sig33(ii) = cf17(xx(ii)) + cf18(xx(ii)) + cf19(xx(ii));
    MS.x(ii) = xx(ii);
    
    MA.u1(ii) = cf21(xx(ii)) + cf22(xx(ii));
    MA.u3(ii) = cf23(xx(ii)) + cf24(xx(ii));
    MA.sig31(ii) = cf25(xx(ii)) + cf26(xx(ii));
    MA.sig33(ii) = cf27(xx(ii)) + cf28(xx(ii)) + cf29(xx(ii));
    MA.x(ii) = xx(ii);
end

end

%----------------------------------------------------------------------------------------------

function [MA MS] = plateModes(w_min, w_max, k_min, k_max)
% Обчислює пластинчаті хвилі як симетричні так і антисиметричні
% головна функція для першого варіанту обчислення, дисперсії хвиль в пластинах

precision = 10;  % поріг, якщо результат функції є меншим то записуємо
dk = 1;    % крок по осі к
[MA] = iterK('acsym', w_min, w_max, k_min, k_max, dk, precision);
precision = 1E-5;  % поріг, якщо результат функції є меншим то записуємо
[MS] = iterK('csym', w_min, w_max, k_min, k_max, dk, precision);
[MA] = cphase(MA);
[MS] = cphase(MS);

end

%----------------------------------------------------------------------------------------------

function [cl ct] = soundVelocity(E, nu, ro)
% Рахує швидкість обємних хвиль за механічними параметрами, модуль Юнга Е, коефіцієнт Пуассона nu
% cl, ct - швидкості, MM.ct, MM.cl  - для виводу

% Коефіцієнти Ляме
mu = E / ( 2 * ( 1 + nu ) );
lambda = (nu*E) / ( (1+nu)*( 1-(2*nu) ) );
cl = sqrt( (lambda + (2*mu) ) / ro );
ct = sqrt( mu / ro );

end

%--------------------------------------------------------------------------%--------------------

function [ VW ] = rapidSearch(fun, w_min, w_max, k, precision)
% Шукає перші точки мод на відрізку w_min .. w_max
% fun - функція для пошуку кореня, fun(k,w).
% w_min - мінімальна частота діапазону пошуку
% w_max - максимальна частота діапазону пошуку
% precision - точність пошуку кореня алгоритмом
% Результати:
% VW - вектор знайдених мод

% пошук по діапазону w_min ... w_max, розбитому на 500 відрізікв
rootsFound = 0; % кількість знайдених коренів
iterw = 0;
stepw = w_max/500;
for ww=w_min:stepw:(w_max-stepw)
    f1 = feval(fun, k, ww);
    f2 = feval(fun, k, ww + stepw);
    if (f1 * f2) < 0
            rootsFound = rootsFound+1;
            iterw(1,rootsFound) = ww;
            iterw(2,rootsFound) = ww + stepw;
    end
end

% уточнення w та виокремлення мод
for iterat = 1:5    % 5 разів уточнююється в 10 раз зменшуючи відрізок, тому 10*10*10*10*10=100 000
    for modeii = 1:rootsFound
        modeLeft = iterw(1, modeii);
        modeRight = iterw(2, modeii);
        stepw = (modeRight - modeLeft) / 10;

        for ww = modeLeft:stepw:(modeRight - stepw)
            f1 = feval(fun, k, ww);
            f2 = feval(fun, k, ww + stepw);
            if (f1 * f2) < 0
                iterw(1,modeii) = ww;
                iterw(2,modeii) = ww + stepw;
                break;
            end
        end

    end
end

VW = zeros(1, rootsFound);
% здійснюємо пошук кореня
for ii=1:rootsFound
    funsearch = @(ww) feval(fun, k, ww);
    [ww] = bisect(funsearch, iterw(1,ii), iterw(2,ii), precision);
    VW(ii) = ww;      % перший рядок w

end

end

%----------------------------------------------------------------------------------------------

function [MM] = iterK(fun, w_min, w_max, k_min, k_max, dk, precision)
% здійснює пошук коренів по осі к для діапазону частот
% Результати
% MM.w - матриця у кожному рядку всі корені для одного значення к
% MM.k - стовпчик значень к

MM.k=0;
MM.w=0;
index = 0;
MM.k = k_min;
for kk = k_min : dk : k_max
    [ VW ] = rapidSearch(fun, w_min, w_max, kk, precision);
    if ~isempty(VW)
        index = index + 1;
        MM.k(index) = kk;
        for ww=1:length(VW)
            MM.w(index, ww) = VW(ww);
        end
    end  
end
MM.k = MM.k';

end

%----------------------------------------------------------------------------------------------

function [MM] = cphase(MM)
% Обчислює фазову швидкість сф = w / k
klen = length(MM.k);
for ii=1:klen
    ww = length(MM.w(ii,:));
    for jj = 1:ww
        MM.cph(ii,jj) = MM.w(ii,jj) / MM.k(ii);
    end
end
end

%----------------------------------------------------------------------------------------------

function [ws] = csym(k,w)
% Обчислення дисперсійних значень симетричних мод
% k - хвильове число
% w - частота
% Результати:
% ws - результат обчислення характеристичного рівняння поширення симетричної моди в пластині.
% Рівняння має дорівнювати нулю, коли мода є.
global d cl ct;
q2=(w/ct)^2-k^2;   % деякі змінні в квадраті
p2=(w/cl)^2-k^2;   % деякі змінні в квадраті
h = d/2;    % 2h=d товщина стінки пластини
if k > (w/ct)
    qq = sqrt(q2)/i; pp = sqrt(p2)/i;
    qq2 = -q2;
elseif ((w/ct)>k) && (k>(w/cl))
    qq = sqrt(q2); pp = sqrt(p2)/i;
    qq2 = q2;
else
    qq = sqrt(q2); pp = sqrt(p2);
    qq2 = q2;
end

part1 = tan(qq*h)/qq;
part2 = (4*k*k*pp*tan(pp*h)) / ((qq2-k*k)^2);
ws = part1+part2;
end

%----------------------------------------------------------------------------------------------

function [ws] = acsym(k,w)
% Обчислення дисперсійних значень антисиметричних мод
% k - хвильове число
% w - частота
% Результати:
% ws - результат обчислення характеристичного рівняння поширення антисиметричної моди в пластині.
% Рівняння має дорівнювати нулю, коли мода є.
global d cl ct;
q2=(w/ct)^2-k*k;   % деякі змінні в квадраті
p2=(w/cl)^2-k*k;   % деякі змінні в квадраті
h = d/2;    % 2h=d товщина стінки пластини

if k > (w/ct)
    qq = sqrt(q2)/i; pp = sqrt(p2)/i;
    qq2 = -q2;
elseif ((w/ct)>k) && (k>(w/cl))
    qq = sqrt(q2); pp = sqrt(p2)/i;
    qq2 = q2;
else
    qq = sqrt(q2); pp = sqrt(p2);
    qq2 = q2;
end

part1 = qq*tan(qq*h);
part2 = (((qq2-k*k)^2)*tan(pp*h)) / (4*k*k*pp);
ws = part1+part2;
end

%----------------------------------------------------------------------------------------------

function [c,yc] = bisect(f,a,b,delta)
%BISECT   The bisection method is used to locate a root.
% Sample calls
%   [c,yc] = bisect('f',a,b,delta)
% Inputs
%   f       name of the function
%   a       left endpoint
%   b       right endpoint
%   delta   convergence tolerance
% Результати:
%   c       solution: the root
%   yc      solution: the function value
%

ya = feval(f,a);
yb = feval(f,b);

if ya > 0
    snga = 1;
elseif ya < 0
    snga = -1;
else
    snga = 0;
end

if yb > 0
    sngb = 1;
elseif yb < 0
    sngb = -1;
else
    sngb = 0;
end

% якщо немає зміни знака то виходимо
if snga*sngb > 0
    yc = 1234512345;
    c = 0;
    return;
end

if ya > 0
    c=b;
    b=a;
    a=c;
    yb = feval(f,b);
end

yc = delta + 10;
% перевірка: якщо через checkNN кроків не зменшується ya або yb то виходимо
checkNN = 5;
checkYA = abs(feval(f,a));
checkYB = abs(feval(f,b));
while abs(yc) > delta
  c  = (a+b)/2;
  yc = feval(f,c);
  
  if yb > 0
      sngb = 1;
  elseif yb < 0
      sngb = -1;
  else
      sngb = 0;
  end
  
  if yc > 0
      sngc = 1;
  elseif yc < 0
      sngc = -1;
  else
      sngc = 0;
  end
  
  if  yc == 0
    a = c;
    b = c;
  elseif  sngb*sngc > 0
    b = c;
    yb = yc;
  else
    a = c;
  end
  % перевірка на сходжуваність
  checkNN = checkNN - 1;
  if checkNN == 0
      if ( checkYA <= abs(feval(f,a)) ) && ( checkYB <= abs(feval(f,b)) )
          return;
      end
  end
  
end
end

%----------------------------------------------------------------------------------------------

