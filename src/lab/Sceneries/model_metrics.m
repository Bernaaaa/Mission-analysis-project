function [rmse, r2] = model_metrics(predicted, actual)
    % Lavora su vettori, ignora NaN
    p = predicted(:);
    a = actual(:);
    valid = ~isnan(p) & ~isnan(a);
    p = p(valid); a = a(valid);
    
    residuals = a - p;
    ss_res    = sum(residuals.^2);
    ss_tot    = sum((a - mean(a)).^2);
    
    rmse = sqrt(mean(residuals.^2));
    r2   = 1 - ss_res / ss_tot;
end