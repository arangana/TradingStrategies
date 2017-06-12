function [p, settings] = ADX(DATE, OPEN, HIGH, LOW, CLOSE, VOL, exposure, equity, settings)
    % Trading System Attributes
    settings.markets     = {'CASH', 'F_AD', 'F_BO', 'F_BP', 'F_C', 'F_CC', 'F_CD', 'F_CL', 'F_CT', 'F_DX', 'F_EC', 'F_ED', 'F_ES', 'F_FC', 'F_FV', 'F_GC', 'F_HG', 'F_HO', 'F_JY', 'F_KC', 'F_LB', 'F_LC', 'F_LN', 'F_MD', 'F_MP', 'F_NG', 'F_NQ', 'F_NR', 'F_O', 'F_OJ', 'F_PA', 'F_PL', 'F_RB', 'F_RU', 'F_S', 'F_SB', 'F_SF', 'F_SI', 'F_SM', 'F_TU', 'F_TY', 'F_US', 'F_W', 'F_XX', 'F_YM'};
    %settings.samplebegin = 20100321;
    %settings.sampleend = 20170320;
    settings.budget = 1000000;
    settings.slippage = 0.05;
    settings.lookback = 200;    
    
    %%% Custom Fields %%%
    settings.num_markets = length(settings.markets);
    if (not(isfield(settings, 'TradeDay')))
        settings.TradeDay = 0;
        % 14 day smoothed True Range
        settings.TR14 = zeros(1,settings.num_markets);
        % 14-day smoothed Positive Directional Movement
        settings.plus_DM14 = zeros(1,settings.num_markets);
        % 14-day smoothed Positive Directional Movement
        settings.minus_DM14 = zeros(1,settings.num_markets);
        % ADX values for each of the selected markets ordered by index
        settings.ADX = NaN(1,settings.num_markets);
    end
    
    %{ 
    ADX Period - 14 is the most commonly used. Set as necessary - 100 is max.
    %}
    ADX_period = 30; %#[7:7:100]#
    
    
    settings.TradeDay = settings.TradeDay + 1;
    disp(['Executing trade for day ', num2str(settings.TradeDay), ' on trade date ', num2str(DATE(settings.lookback))]);

    % Compute ADX for each market
    for market = 1:settings.num_markets
        endDate = settings.lookback;
        startDate = endDate - endDate + 1;
        settings = calculate_ADX(market, HIGH(startDate:endDate,market), LOW(startDate:endDate,market), CLOSE(startDate:endDate,market), settings, ADX_period);
    end
    
    % Execute trades based on trading strategy
    p = execute_trade(settings.ADX, CLOSE, ADX_period);
     
end

function weights = execute_trade(market_ADXs, CLOSE, ADX_period)

    CURRENT_DAY = length(CLOSE(:,1));
    num_markets = length(market_ADXs);
    ADX_sum = nansum(market_ADXs);
    % Proportion of capital for market exposure for each equity
    prop = zeros(1,num_markets);
    
    % 14 day moving average of closing prices
    endDay = (ADX_period * 2) - 1;
    startDay = endDay - ADX_period + 1;
    avg = sum(CLOSE(startDay:endDay,:)) / ADX_period;
    
    for market = 1:num_markets
        % There is a strong trend
        if (market_ADXs(market) > 25.0)
            % Strong trend indicating rising market - BUY
            if (CLOSE(CURRENT_DAY,market) > avg(market))
                prop(market) = 1.0 * market_ADXs(market) / ADX_sum;
            % Strong trend indicating falling market - SELL
            elseif (CLOSE(CURRENT_DAY,market) < avg(market))
                prop(market) = -1.0 * market_ADXs(market) / ADX_sum;
            else
                prop(market) = 0.0;
            end
        % No strong trend, no market exposure
        else
            prop(market) = 0.0;
        end  
    end

    weights = prop;
    
end

% Details of ADX calculations for a given market
function ret = calculate_ADX(market, HIGH, LOW, CLOSE, settings, ADX_period)

    % Components required to calculate ADX - True Range, Directional Movements and Directional Indices
    tradeDayVals.TR = 0.0;
    tradeDayVals.plus_DM = 0.0;
    tradeDayVals.minus_DM = 0.0;
    tradeDayVals.plus_DI = 0.0;
    tradeDayVals.minus_DI = 0.0;
    tradeDayVals.DX = 0.0;
    
    % IPO has not happened long enough ago to calculate ADX
    if (isnan(CLOSE(1)))
        settings.ADX(market) = NaN;
    
    % First calculation of ADX - i.e. no previous values of ADX available
    elseif (isnan(settings.ADX(market)))
        for trade_day = 2:(ADX_period + 1)
            % Add up / accumulate true range and directional movement values for first 14 days  
            settings.TR14(market)= settings.TR14(market) + calculate_TR(CLOSE(trade_day-1), HIGH(trade_day), LOW(trade_day));
            settings.plus_DM14(market) = settings.plus_DM14(market) + calculate_DM_plus(HIGH(trade_day), LOW(trade_day), HIGH(trade_day-1), LOW(trade_day-1));
            settings.minus_DM14(market) = settings.minus_DM14(market) + calculate_DM_minus(HIGH(trade_day), LOW(trade_day), HIGH(trade_day-1), LOW(trade_day-1));
        end
        
        % First calculation of directional indices
        tradeDayVals.plus_DI = 100.0 * settings.plus_DM14(market) / settings.TR14(market);
        tradeDayVals.minus_DI = 100.0 * settings.minus_DM14(market) / settings.TR14(market);
        settings.ADX(market) = calculate_DX(tradeDayVals.plus_DI, tradeDayVals.minus_DI);
        
        for trade_day = (ADX_period + 2):(ADX_period * 2)

            % True range and directional movements for the current trade day
            tradeDayVals.TR = calculate_TR(CLOSE(trade_day-1), HIGH(trade_day), LOW(trade_day));
            tradeDayVals.plus_DM = calculate_DM_plus(HIGH(trade_day), LOW(trade_day), HIGH(trade_day-1), LOW(trade_day-1));
            tradeDayVals.minus_DM = calculate_DM_minus(HIGH(trade_day), LOW(trade_day), HIGH(trade_day-1), LOW(trade_day-1));

            % Wilder smoothing techniques to calculate 14 day smoothed true range and directional movements
            settings.TR14(market) = wilder_smooth(settings.TR14(market), tradeDayVals.TR, ADX_period);
            settings.plus_DM14(market) = wilder_smooth(settings.plus_DM14(market), tradeDayVals.plus_DM, ADX_period);
            settings.minus_DM14(market) = wilder_smooth(settings.minus_DM14(market), tradeDayVals.minus_DM, ADX_period); 

            % Directional indices for the current trade day
            tradeDayVals.plus_DI = 100.0 * settings.plus_DM14(market) / settings.TR14(market);
            tradeDayVals.minus_DI = 100.0 * settings.minus_DM14(market) / settings.TR14(market);
            tradeDayVals.DX = calculate_DX(tradeDayVals.plus_DI, tradeDayVals.minus_DI);

            % Add up / accumulate DX values for next 14 days 
            settings.ADX(market) = settings.ADX(market) + tradeDayVals.DX;
            %disp(sprintf('ADX vals accum: (%f, %f, %f, %f, %f)', settings.ADX));
        end
        % First calculation of ADX - average of first 14 DX values
        settings.ADX(market) = settings.ADX(market) / ADX_period;
        
    % Subsequent calculations of ADX
    else
        CURRENT_DAY = length(HIGH);
        PREVIOUS_DAY = CURRENT_DAY - 1;
        
        % True Range and Directional Movement Calculations for current trading day
        tradeDayVals.TR = calculate_TR(CLOSE(PREVIOUS_DAY), HIGH(CURRENT_DAY), LOW(CURRENT_DAY));
        tradeDayVals.plus_DM = calculate_DM_plus(HIGH(CURRENT_DAY), LOW(CURRENT_DAY), HIGH(PREVIOUS_DAY), LOW(PREVIOUS_DAY));
        tradeDayVals.minus_DM = calculate_DM_minus(HIGH(CURRENT_DAY), LOW(CURRENT_DAY), HIGH(PREVIOUS_DAY), LOW(PREVIOUS_DAY));

        % Wilder smoothing techniques to calculate 14 day smoothed true range and directional movements
        settings.TR14(market) = wilder_smooth(settings.TR14(market), tradeDayVals.TR, ADX_period);
        settings.plus_DM14(market) = wilder_smooth(settings.plus_DM14(market), tradeDayVals.plus_DM, ADX_period);
        settings.minus_DM14(market) = wilder_smooth(settings.minus_DM14(market), tradeDayVals.minus_DM, ADX_period); 
        
        % Directional indices for the current trade day
        tradeDayVals.plus_DI = 100.0 * settings.plus_DM14(market) / settings.TR14(market);
        tradeDayVals.minus_DI = 100.0 * settings.minus_DM14(market) / settings.TR14(market);
        tradeDayVals.DX = calculate_DX(tradeDayVals.plus_DI, tradeDayVals.minus_DI);
        
        % Final ADX value for current trade day
        settings.ADX(market) = ((settings.ADX(market) * (ADX_period-1)) + (tradeDayVals.DX)) / ADX_period;
        
    end
    
    ret = settings;
    
end


function ret = wilder_smooth(prev_val, curr_val, ADX_period)

    ret = prev_val - (prev_val / ADX_period) + curr_val;

end

function ret = calculate_TR(prevClose, currHigh, currLow)

    criterion_1 = currHigh - currLow;
    criterion_2 = abs(currHigh - prevClose);
    criterion_3 = abs(prevClose - currLow);
    ret = max([criterion_1, criterion_2, criterion_3]);

end

function ret = calculate_DM_plus(currHigh, currLow, prevHigh, prevLow)

    currHigh_minus_prevHigh = currHigh - prevHigh;
    if (currHigh_minus_prevHigh <= 0)
        ret = 0.0;
        return;
    end
    prevLow_minus_currLow = prevLow - currLow;
    if(currHigh_minus_prevHigh > prevLow_minus_currLow)
        ret = currHigh_minus_prevHigh;
        return;
    end
    ret =  0.0;
 
end

function ret = calculate_DM_minus(currHigh, currLow, prevHigh, prevLow)

    prevLow_minus_currLow = prevLow - currLow;
    if (prevLow_minus_currLow <= 0)
        ret = 0.0;
        return;
    end
    currHigh_minus_prevHigh = currHigh - prevHigh;
    if(prevLow_minus_currLow > currHigh_minus_prevHigh)
        ret =  prevLow_minus_currLow;
        return;
    end
    ret = 0.0;

end

function ret = calculate_DX(plus_DI, minus_DI)

    ret = 100.0 * (abs(plus_DI - minus_DI)) / (plus_DI + minus_DI);

end