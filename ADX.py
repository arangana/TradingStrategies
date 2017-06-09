### Simple trading strategy based on ADX calculations. ###

import numpy as np


''' This system uses Average Directional Index (ADX) to allocate capital across desired equities. '''
def myTradingSystem(DATE, OPEN, HIGH, LOW, CLOSE, VOL, exposure, equity, settings):

    settings['TradeDay'] = settings['TradeDay'] + 1
    print 'Executing trade for day ', settings['TradeDay'], ' on trade date ', DATE[settings['lookback']-1]

    # Compute ADX for each market
    num_markets = settings['num_markets']
    for market in range(num_markets):
        calculate_ADX(market, HIGH[:,market:market+1], LOW[:,market:market+1], CLOSE[:,market:market+1], settings)

    # Execute trades based on trading strategy
    weights = execute_trade(settings['ADX'], CLOSE)

    return weights, settings


''' Trading system settings. '''
def mySettings():

    settings= {}

    settings['beginInSample'] = '20160321'
    settings['endInSample'] = '20170320'
    settings['markets']  = ['AAPL', 'BRKB', 'F_NG', 'F_GC']
    settings['lookback'] = 28
    settings['budget'] = 10**6
    settings['slippage'] = 0.05    

    ### Custom Fields ###
    num_markets = len(settings['markets'])
    settings['num_markets'] = num_markets
    settings['TradeDay'] = 0
    # 14 day smoothed True Range
    settings['TR14'] = [0.0] * settings['num_markets']
    # 14-day smoothed Positive Directional Movement
    settings['+DM14'] = [0.0] * settings['num_markets']
    # 14-day smoothed Positive Directional Movement
    settings['-DM14'] = [0.0] * settings['num_markets']
    # ADX values for each of the selected markets ordered by index
    settings['ADX'] = [np.nan] * settings['num_markets']

    return settings


'''
Simple trading strategy based on ADX. If ADX value is above 25, we
know that there is a strong trend. We then look at the average price
over the last 14 trading days and check if there is an upward or 
downward trend. Based on this we take on either a long or short position accordingly.
'''
def execute_trade(market_ADXs, CLOSE):
    CURRENT_DAY = CLOSE.shape[0] - 1
    num_markets = len(market_ADXs)
    ADX_sum = np.nansum(market_ADXs)
    # Proportion of capital for market exposure for each equity
    prop = [0.0] * num_markets

    # Prior 14 day average of closing prices
    avg = np.nansum(CLOSE[13:27,:],axis=0) / 14.0

    for market in range(num_markets):
        # There is a strong trend
        if (market_ADXs[market] > 25.0):
            # Strong trend indicating rising market - BUY
            if (CLOSE[CURRENT_DAY][market] > avg[market]):
                prop[market] = 1.0 * market_ADXs[market] / ADX_sum
            # Strong trend indicating falling market - SELL
            elif (CLOSE[CURRENT_DAY][market] < avg[market]):
                prop[market] =  -1.0 * market_ADXs[market] / ADX_sum
            else:
                prop[market] = 0.0
        # No strong trend, no market exposure
        else:
            prop[market] = 0.0
    weights = np.ndarray((num_markets,), buffer=np.array(prop), dtype=float)
    return weights


# Details of ADX calculations for a given market
def calculate_ADX(market, HIGH, LOW, CLOSE, settings):
    # Components required to calculate ADX - True Range, Directional Movements and Directional Indices
    tradeDayVals = {'TR':0.0, '+DM':0.0, '-DM':0.0, '+DI':0.0, '-DI':0.0, 'DX':0.0}

    # IPO has not happened long enough ago to calculate ADX
    if(np.isnan(CLOSE[0][0])):
        settings['ADX'][market] = np.nan
        
    # First calculation of ADX - i.e. no previous values of ADX available
    elif(np.isnan(settings['ADX'][market])):
        for trade_day in range(1, 15):
            # Add up / accumulate true range and directional movement values for first 14 days         
            settings['TR14'][market] = settings['TR14'][market] + calculate_TR(CLOSE[trade_day-1][0], HIGH[trade_day][0], LOW[trade_day][0])
            settings['+DM14'][market] = settings['+DM14'][market] + calculate_DM_plus(HIGH[trade_day][0], LOW[trade_day][0], HIGH[trade_day-1][0], LOW[trade_day-1][0])
            settings['-DM14'][market] = settings['-DM14'][market] + calculate_DM_minus(HIGH[trade_day][0], LOW[trade_day][0], HIGH[trade_day-1][0], LOW[trade_day-1][0])

        # First calculation of directional indices
        tradeDayVals['+DI'] = 100.0 * settings['+DM14'][market] / settings['TR14'][market]
        tradeDayVals['-DI'] = 100.0 * settings['-DM14'][market] / settings['TR14'][market]
        settings['ADX'][market] = calculate_DX(tradeDayVals['+DI'], tradeDayVals['-DI'])

        for trade_day in range(15, 28):
            # True range and directional movements for the current trade day
            tradeDayVals['TR'] = calculate_TR(CLOSE[trade_day-1][0], HIGH[trade_day][0], LOW[trade_day][0])
            tradeDayVals['+DM'] = calculate_DM_plus(HIGH[trade_day][0], LOW[trade_day][0], HIGH[trade_day-1][0], LOW[trade_day-1][0])
            tradeDayVals['-DM'] = calculate_DM_minus(HIGH[trade_day][0], LOW[trade_day][0], HIGH[trade_day-1][0], LOW[trade_day-1][0])

            # Wilder smoothing techniques to calculate 14 day smoothed true range and directional movements
            settings['TR14'][market] = wilder_smooth(settings['TR14'][market], tradeDayVals['TR'])
            settings['+DM14'][market] = wilder_smooth(settings['+DM14'][market], tradeDayVals['+DM'])
            settings['-DM14'][market] = wilder_smooth(settings['-DM14'][market], tradeDayVals['-DM'])

            # Directional indices for the current trade day
            tradeDayVals['+DI'] = 100.0 * settings['+DM14'][market] / settings['TR14'][market]
            tradeDayVals['-DI'] = 100.0 * settings['-DM14'][market] / settings['TR14'][market]
            tradeDayVals['DX'] = calculate_DX(tradeDayVals['+DI'], tradeDayVals['-DI'])

            # Add up / accumulate DX values for next 14 days 
            settings['ADX'][market] = settings['ADX'][market] + tradeDayVals['DX']

        # First calculation of ADX - average of first 14 DX values
        settings['ADX'][market] = settings['ADX'][market] / 14.0

    # Subsequent calculations of ADX
    else:
        CURRENT_DAY = HIGH.shape[0] - 1
        PREVIOUS_DAY = CURRENT_DAY - 1

        # True Range and Directional Movement Calculations for current trading day
        tradeDayVals['TR'] = calculate_TR(CLOSE[PREVIOUS_DAY][0], HIGH[CURRENT_DAY][0], LOW[CURRENT_DAY][0])
        tradeDayVals['+DM'] = calculate_DM_plus(HIGH[CURRENT_DAY][0], LOW[CURRENT_DAY][0], HIGH[PREVIOUS_DAY][0], LOW[PREVIOUS_DAY][0])
        tradeDayVals['-DM'] = calculate_DM_minus(HIGH[CURRENT_DAY][0], LOW[CURRENT_DAY][0], HIGH[PREVIOUS_DAY][0], LOW[PREVIOUS_DAY][0])

        # Wilder smoothing techniques to calculate 14 day smoothed true range and directional movements
        settings['TR14'][market] = wilder_smooth(settings['TR14'][market], tradeDayVals['TR'])
        settings['+DM14'][market] = wilder_smooth(settings['+DM14'][market], tradeDayVals['+DM'])
        settings['-DM14'][market] = wilder_smooth(settings['-DM14'][market], tradeDayVals['-DM'])

        # Directional indices for the current trade day
        tradeDayVals['+DI'] = 100.0 * settings['+DM14'][market] / settings['TR14'][market]
        tradeDayVals['-DI'] = 100.0 * settings['-DM14'][market] / settings['TR14'][market]
        tradeDayVals['DX'] = calculate_DX(tradeDayVals['+DI'], tradeDayVals['-DI'])

        # Final ADX value for current trade day
        settings['ADX'][market] = ((settings['ADX'][market] * 13.0) + (tradeDayVals['DX'])) / 14.0

    return


def wilder_smooth(prev_val, curr_val):
    return prev_val - (prev_val / 14.0) + curr_val


def calculate_TR(prevClose, currHigh, currLow):
    criterion_1 = currHigh - currLow
    criterion_2 = abs(currHigh - prevClose)
    criterion_3 = abs(prevClose - currLow)
    return max(criterion_1, criterion_2, criterion_3)


def calculate_DM_plus(currHigh, currLow, prevHigh, prevLow):
    currHigh_minus_prevHigh = currHigh - prevHigh
    if (currHigh_minus_prevHigh <= 0):
        return 0.0
    prevLow_minus_currLow = prevLow - currLow
    if(currHigh_minus_prevHigh > prevLow_minus_currLow):
        return currHigh_minus_prevHigh
    return 0.0


def calculate_DM_minus(currHigh, currLow, prevHigh, prevLow):
    prevLow_minus_currLow = prevLow - currLow
    if (prevLow_minus_currLow <= 0):
        return 0.0
    currHigh_minus_prevHigh = currHigh - prevHigh
    if(prevLow_minus_currLow > currHigh_minus_prevHigh):
        return prevLow_minus_currLow
    return 0.0
    
    
def calculate_DX(plus_DI, minus_DI):
    return 100.0 * (abs(plus_DI - minus_DI)) / (plus_DI + minus_DI)


if __name__ == '__main__':
    import quantiacsToolbox
    results = quantiacsToolbox.runts(__file__)
