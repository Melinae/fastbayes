{-# LANGUAGE UnicodeSyntax #-}
{-# LANGUAGE BangPatterns #-}

module Statistics.FastBayes.Linear 
  (Fit, marginalLikelihood)
  where

import qualified Data.Vector.Storable as V
import Numeric.LinearAlgebra

α0 :: Double
α0 = 1.0

square :: Double -> Double
square x = x * x

normSq x = x <.> x

data Fit = Fit
  { design                 :: Matrix Double
  , response               :: Vector Double
  , priorPrecision         :: Double
  , noisePrecision         :: Double
  , numEffectiveParameters :: Double
  , logEvidence            :: Double
  , mapWeights             :: Vector Double
  , hessian                :: Matrix Double
  }
  deriving Show


marginalLikelihood :: 
  ([(Double, Double)] → (Double, Double))  -- limit of the (α,β) sequence
  → Matrix Double                           -- design matrix (features in columns)
  → Vector Double                           -- response vector
  → Fit
marginalLikelihood lim x y = Fit x y α β γ logEv m h
  where
  n = rows x
  p = cols x
  α0 = 1.0
  β0 = 1.0

  -- A vector of the eigenvalues of xtx
  eig = V.map square $ singularValues x

  xtx = trans x <> x
  xty = trans x <> y
  hessian α β = diag (V.replicate p α) + scale β xtx

  m = scale β $ h <\> xty
  
  go :: Double → Double → [(Double, Double)]
  go α0 β0 = (α0, β0) : go α β
    where
    h = hessian α0 β0
    m = scale β0 $ h <\> xty 
    γ = V.sum $ V.map (\x → x / (α0 + x)) eig
    α = γ / (m <.> m)
    β = recip $ (normSq $ y - x <> m) / (fromIntegral n - γ)
  
  γ = V.sum $ V.map (\x → x / (α + x)) eig

  h = hessian α β 
  logEv = 0.5 * 
    ( fromIntegral p * log α 
    + fromIntegral n * log β 
    - (β * normSq (y - x <> m) + α * (m <.> m))
    - logDetH 
    - fromIntegral n * log (2*pi)
    )
    where
    (_,(logDetH, _)) = invlndet h

  (α, β) = lim $ go α0 β0


x1 = gaussianSample 3 1000 (V.fromList $ replicate 100 0) (ident 100)
x2 = trans $ (100 >< 1000) [1.0 ..]

y = V.fromList [1.0 .. 1000.0]

--x = trans $ (200 >< 5000) [1.0 ..]
--y = V.fromList [1.0 .. 5000.0]

fit1 = marginalLikelihood (!! 100) x1 y
fit2 = marginalLikelihood (!! 100) x2 y
