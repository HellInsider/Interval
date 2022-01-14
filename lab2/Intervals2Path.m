function [P] = Intervals2Path(S)
% 
% Функция Intervals2Path по непустой матрице 
% граничных интервалов S с конечными элементами
% строит замкнутый обход P границы политопа.

    bp = [S(1,1) S(1,2)]; % начало обхода
    P = bp;
    bs = bp; % начало прямолинейного участка границы

    while size(S,1)>0

       for k=1:size(S,1); % ищем строку i, для которой граничный интервал 
                          % S(i,:) начинается в bs
%         if isequal(bs,[S(k,1) S(k,2)])
          if max(abs(bs-[S(k,1) S(k,2)]))<1.e-8
             i=k;
             break;
          end
       end

       es = [S(i,3) S(i,4)]; % конец прямолинейного участка границы
%      if ~isequal(es,bs)
       if max(abs(bs-es))>1.e-8
          P = [P; es];
%         if isequal(es,bp)
          if max(abs(bp-es))<1.e-8
             return  % построен замкнутый путь
          end
          bs = es; % начинаем новый прямолинейный участок границы 
                   % из конца предыдущего
       end
       S(i,:)=[];

    end
    
end