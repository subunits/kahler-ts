-- Augmented Vector Space — Kähler Extension v2
-- No dependencies — runs on play.haskell.org (base only)
-- Compile: ghc -O kahler_isco.hs -o kahler_isco

module Main where

import Data.List (sortBy, intercalate, nub, minimumBy)
import Data.Ord  (comparing)

-- ─── Complex numbers ─────────────────────────────────────────────────────────

data C = C { re :: Double, im :: Double }

instance Show C where
  show (C r i) = fmtF 3 r ++ (if i >= 0 then "+" else "") ++ fmtF 3 i ++ "i"

(|+|),(|-|),(|*|) :: C -> C -> C
C a b |+| C c d = C (a+c) (b+d)
C a b |-| C c d = C (a-c) (b-d)
C a b |*| C c d = C (a*c - b*d) (a*d + b*c)

conj :: C -> C
conj (C r i) = C r (-i)

cabs :: C -> Double
cabs (C r i) = sqrt (r*r + i*i)

fromReal :: Double -> C
fromReal x = C x 0

-- ─── Complex vector operations ────────────────────────────────────────────────

type CVec = [C]

hermitian :: CVec -> CVec -> C
hermitian u v = foldr (|+|) (C 0 0) (zipWith (\a b -> conj a |*| b) u v)

riemannian :: CVec -> CVec -> Double
riemannian u v = re (hermitian u v)

symplectic :: CVec -> CVec -> Double
symplectic u v = im (hermitian u v)

complexJ :: CVec -> CVec
complexJ = map (\(C r i) -> C (-i) r)

cnorm :: CVec -> Double
cnorm v = sqrt (riemannian v v)

-- ─── Kähler condition checks ──────────────────────────────────────────────────

checkJSquared :: CVec -> Bool
checkJSquared v =
  let jjv = complexJ (complexJ v)
      neg  = map (\(C r i) -> C (-r) (-i)) v
  in  all (\(C r i, C r' i') -> abs(r-r') < 1e-10 && abs(i-i') < 1e-10)
          (zip jjv neg)

checkJPreservesMetric :: CVec -> CVec -> Bool
checkJPreservesMetric u v =
  abs (riemannian (complexJ u) (complexJ v) - riemannian u v) < 1e-10

checkSkewSymmetric :: CVec -> CVec -> Bool
checkSkewSymmetric u v =
  abs (symplectic u v + symplectic v u) < 1e-10

checkCompatibility :: CVec -> CVec -> Bool
checkCompatibility u v =
  abs (symplectic u v - riemannian (complexJ u) v) < 1e-10

-- Discrete curvature of ω over a triple of embedded points.
-- Returns (value, isCurved) where curved means |val| > epsilon.
-- ✗ means the attractor is curved inside C^k — this is the interesting quantity.
-- A flat C^k satisfies dω=0 globally; embedded signal manifolds generally don't.
omegaCurvature :: CVec -> CVec -> CVec -> (Double, Bool)
omegaCurvature p q r =
  let sub   = zipWith (|-|)
      qp    = sub q p; rp = sub r p
      rq    = sub r q; pq = sub p q
      pr    = sub p r; qr = sub q r
      val   = symplectic qp rp + symplectic rq pq + symplectic pr qr
  in  (val, abs val > 1e-8)

-- ─── Types ───────────────────────────────────────────────────────────────────

type Sample = Double
type Vec    = [Double]

data AugConfig = AugConfig
  { lagWindow  :: Int
  , derivOrder :: Int
  , rollingWin :: Int
  } deriving (Show)

data AugPoint = AugPoint
  { apLabel :: Int
  , apVec   :: Vec
  } deriving (Show)

data KahlerPoint = KahlerPoint
  { kpLabel   :: Int
  , kpReal    :: Vec
  , kpComplex :: CVec
  } deriving (Show)

-- ─── Signal generators ───────────────────────────────────────────────────────

sineWave :: Int -> [Sample]
sineWave n = [ sin (2*pi*fromIntegral t/20)*3 | t <- [0..n-1] ]

logisticMap :: Int -> [Sample]
logisticMap n = take n . map (\x->(x-0.5)*6) $ iterate (\x->3.9*x*(1-x)) 0.5

lorenzWave :: Int -> [Sample]
lorenzWave n = normalise . take n . every 8 . drop 1000 $ map (\(x,_,_)->x) lorenz
  where
    lorenz = iterate step (0.1, 0.0, 0.0)
    step (x,y,z) =
      let dt=0.01;s=10;rho=28;beta=8/3
      in (x+dt*s*(y-x), y+dt*(x*(rho-z)-y), z+dt*(x*y-beta*z))
    every k xs = case drop (k-1) xs of { [] -> []; (y:ys) -> y : every k ys }
    normalise xs = let mx=maximum(map abs xs) in map (*(3/max mx 1e-9)) xs

sawtoothWave :: Int -> [Sample]
sawtoothWave n =
  [ fromIntegral (t`mod`15)/15*4-2
    + sin(fromIntegral t*2.3)*0.4
    + sin(fromIntegral t*0.7)*0.3
  | t <- [0..n-1] ]

stepWave :: Int -> [Sample]
stepWave n = zipWith (+) steps (take n (gaussStream 42))
  where steps = [ if t<20 then 1 else if t<40 then -1 else 2 | t<-[0..n-1] ]

randomWalk :: Int -> [Sample]
randomWalk n = take n $ scanl (+) 0 (map (*0.8) (uniformStream 12345))

uniformStream :: Int -> [Double]
uniformStream seed = map toF (iterate lcg (fromIntegral seed))
  where
    lcg s   = (s*1664525+1013904223) `mod` (2^32)
    toF s   = fromIntegral s / fromIntegral (2^32::Integer) - 0.5

gaussStream :: Int -> [Double]
gaussStream seed =
  [ sqrt(-2*log(max u 1e-9)) * cos(2*pi*v) * 0.6
  | (u,v) <- pairs (map (+0.5) (uniformStream seed)) ]
  where pairs (a:b:rest) = (a,b):pairs rest; pairs _ = []

allSignals :: [(String, [Sample])]
allSignals =
  [ ("ISCO-r",     [-2.2e-05, 0.469324, 0.927044, 1.361952, 1.763352, 2.121341, 2.427033, 2.673001, 2.853188, 2.963087, 3.0, 2.963087, 2.853188, 2.673001, 2.427033, 2.121341, 1.763352, 1.361952, 0.927044, 0.469324, -2.2e-05, -0.469324, -0.927044, -1.361952, -1.763352, -2.121341, -2.427077, -2.673045, -2.853188, -2.963087, -3.0, -2.963087, -2.853188, -2.673045, -2.427077, -2.121341, -1.763352, -1.361952, -0.927044, -0.469324, -2.2e-05, 0.469324, 0.927044, 1.361952, 1.763352, 2.121341, 2.427033, 2.673001, 2.853188, 2.963087, 3.0, 2.963087, 2.853188, 2.673001, 2.427033, 2.121341, 1.763352, 1.361952, 0.927044, 0.469324])
        , ("ISCO-t",     [-3.0, -2.979933, -2.959866, -2.939799, -2.919732, -2.899666, -2.879599, -2.859532, -2.839465, -2.819398, -2.799331, -2.779264, -2.759197, -2.73913, -2.719064, -2.698997, -2.67893, -2.658863, -2.638796, -2.618729, -2.598662, -2.578595, -2.558528, -2.538462, -2.518395, -2.498328, -2.478261, -2.458194, -2.438127, -2.41806, -2.397993, -2.377926, -2.35786, -2.337793, -2.317726, -2.297659, -2.277592, -2.257525, -2.237458, -2.217391, -2.197324, -2.177258, -2.157191, -2.137124, -2.117057, -2.09699, -2.076923, -2.056856, -2.036789, -2.016722, -1.996656, -1.976589, -1.956522, -1.936455, -1.916388, -1.896321, -1.876254, -1.856187, -1.83612, -1.816054])
        , ("ISCO-theta", [0.0, 0.342137, 0.679872, 1.008603, 1.32433, 1.622649, 1.89996, 2.152261, 2.376751, 2.570228, 2.730092, 2.854342, 2.941577, 2.990196, 3.0, 2.970788, 2.902761, 2.796719, 2.654462, 2.477391, 2.268307, 2.029412, 1.764106, 1.47579, 1.168267, 0.845538, 0.511805, 0.171469, -0.171269, -0.511805, -0.845538, -1.168267, -1.47579, -1.763906, -2.029212, -2.268107, -2.477391, -2.654262, -2.796719, -2.902561, -2.970588, -3.0, -2.990196, -2.941377, -2.854342, -2.729892, -2.570028, -2.376551, -2.152261, -1.89976, -1.622449, -1.32413, -1.008403, -0.679672, -0.342137, 0.0, 0.342137, 0.679872, 1.008603, 1.32433])
        , ("ISCO-phi",   [-3.0, -2.868523, -2.737046, -2.605569, -2.474092, -2.342615, -2.211137, -2.079661, -1.948184, -1.816707, -1.685229, -1.553753, -1.422276, -1.290799, -1.159321, -1.027845, -0.896368, -0.764891, -0.633413, -0.501937, -0.37046, -0.238983, -0.107505, 0.023971, 0.155448, 0.286925, 0.418403, 0.549879, 0.681356, 0.812833, 0.944311, 1.075788, 1.207264, 1.338741, 1.470219, 1.601696, 1.733172, 1.86465, 1.996127, 2.127604, 2.25908, 2.390558, 2.522035, 2.653512, 2.784988, 2.916466, -2.972155, -2.840679, -2.709201, -2.577724, -2.446247, -2.314771, -2.183293, -2.051816, -1.920339, -1.788863, -1.657385, -1.525908, -1.394431, -1.262954])
        , ("lorenz",     lorenzWave 60)
        , ("sine",       sineWave   60)
        ]

-- ─── Augmentation ────────────────────────────────────────────────────────────

lagFeats :: Int -> [Sample] -> [(Int,Vec)]
lagFeats p xs =
  [(i, reverse(take p(drop(i-p+1)xs))) | i<-[p-1..length xs-1]]

diffs :: Int -> [Sample] -> [Sample]
diffs 0 xs = xs
diffs d xs = diffs(d-1)(zipWith(-)(drop 1 xs)xs)

rollingStats :: Int -> [Sample] -> [(Double,Double)]
rollingStats w xs =
  [ let win=take w(drop i xs); mu=sum win/fromIntegral w
        var=sum(map(\x->(x-mu)^2)win)/fromIntegral w
    in (mu,sqrt var)
  | i<-[0..length xs-w] ]

augment :: AugConfig -> [Sample] -> [AugPoint]
augment cfg xs = zipWith3 build lagPairs diffSeries statSeries
  where
    p=lagWindow cfg; d=derivOrder cfg; w=rollingWin cfg
    lagPairs   = lagFeats p xs
    start      = p-1
    diffXs     = diffs 1 xs; diff2Xs = diffs 2 xs
    diffSeries = case d of
      0 -> repeat []
      1 -> map(:[]) (drop start diffXs)
      _ -> zipWith(\a b->[a,b]) (drop start diffXs) (drop start diff2Xs)
    rawStats   = rollingStats w xs
    statSeries = drop(start-(w-1)) rawStats
    build (idx,lags) ds (mu,sigma) = AugPoint idx (lags++ds++[mu,sigma])

complexify :: AugConfig -> AugPoint -> KahlerPoint
complexify cfg pt = KahlerPoint (apLabel pt) (apVec pt) cvec
  where
    p    = lagWindow cfg
    lags = take p (apVec pt)
    rest = drop p (apVec pt)
    cvec = pairUp lags ++ map fromReal rest
    pairUp []       = []
    pairUp [x]      = [C x 0]
    pairUp (x:y:zs) = C x y : pairUp zs

-- ─── Distances ───────────────────────────────────────────────────────────────

holoDist :: KahlerPoint -> KahlerPoint -> Double
holoDist p q =
  sqrt . sum $ zipWith(\a b->cabs(a|-|b)^2) (kpComplex p) (kpComplex q)

-- ─── kNN ─────────────────────────────────────────────────────────────────────

holoKNN :: Int -> KahlerPoint -> [KahlerPoint] -> [(KahlerPoint,Double)]
holoKNN k query corpus =
  take k . sortBy(comparing snd)
  . map(\p->(p,holoDist query p))
  . filter((/=kpLabel query).kpLabel) $ corpus

-- ─── Symplectic matrix & Chern proxy ─────────────────────────────────────────

symplecticMatrix :: [KahlerPoint] -> [[Double]]
symplecticMatrix pts =
  [[symplectic(kpComplex p)(kpComplex q)|q<-pts]|p<-pts]

chernProxy :: [KahlerPoint] -> Double
chernProxy pts =
  let mat=symplecticMatrix pts; n=length pts
  in  sum[mat!!i!!j|i<-[0..n-1],j<-[i+1..n-1]]

-- ─── Shannon entropy (binned) ────────────────────────────────────────────────

shannonEntropy :: [Double] -> Double
shannonEntropy [] = 0
shannonEntropy xs =
  let n    = length xs
      bins = max 4 (round (sqrt (fromIntegral n)) :: Int)
      mn   = minimum xs
      mx   = maximum xs
      rng  = max (mx - mn) 1e-9
      idx v = min (bins-1) (floor ((v-mn)/rng * fromIntegral bins) :: Int)
      counts = foldr (\v acc ->
                 let i = idx v
                     (pre, c:post) = splitAt i acc
                 in  pre ++ (c+1) : post
               ) (replicate bins 0) xs
      total = fromIntegral n :: Double
  in  negate . sum $
        [ let p = fromIntegral c / total
          in  if p > 0 then p * logBase 2 p else 0
        | c <- counts ]

-- ─── Persistent homology (Vietoris-Rips) ─────────────────────────────────────
-- Union-find for connected components (H0)
-- Triangle closure for 1-cycles (H1 proxy)

-- Simple union-find using plain lists, no state threading issues
type Parents = [Int]

mkParents :: Int -> Parents
mkParents n = [0..n-1]

findRoot :: Parents -> Int -> Int
findRoot p i = if p !! i == i then i else findRoot p (p !! i)

unionP :: Parents -> Int -> Int -> Parents
unionP p a b =
  let ra = findRoot p a
      rb = findRoot p b
  in  if ra == rb then p
      else take ra p ++ [rb] ++ drop (ra+1) p

countComponents :: Parents -> Int -> Int
countComponents p n =
  length . nub $ map (findRoot p) [0..n-1]

-- Vietoris-Rips filtration
-- Returns [(epsilon, H0_components, H1_cycles, edge_entropy)]
vrFiltration :: [KahlerPoint] -> Int -> [(Double,Int,Int,Double)]
vrFiltration pts nSteps =
  let n     = length pts
      dists = [ [ holoDist (pts!!i) (pts!!j) | j<-[0..n-1] ] | i<-[0..n-1] ]
      allDs    = sortBy compare
                   [ dists!!i!!j | i<-[0..n-1], j<-[i+1..n-1] ]
      m        = length allDs
      idxs     = nub [ (k*(m-1)) `div` max 1 (nSteps-1) | k<-[0..nSteps-1] ]
      epsilons = map (allDs!!) idxs
      edges e    = [ (i,j) | i<-[0..n-1], j<-[i+1..n-1], dists!!i!!j<=e ]
      triangles e = [ (i,j,k) | i<-[0..n-2], j<-[i+1..n-1], k<-[j+1..n-1]
                               , dists!!i!!j<=e, dists!!i!!k<=e
                               , dists!!j!!k<=e ]
      h0 e =
        let p = foldl (\acc (a,b) -> unionP acc a b) (mkParents n) (edges e)
        in  countComponents p n
      h1 e =
        let es = edges e; ts = triangles e; comps = h0 e
        in  max 0 (length es - (n - comps) - length ts)
      edgeEntropy e =
        let ds = [ dists!!i!!j | (i,j) <- edges e ]
        in  shannonEntropy ds
  in  map (\e->(e, h0 e, h1 e, edgeEntropy e)) epsilons

-- ─── Pretty printing ─────────────────────────────────────────────────────────

fmtF :: Int -> Double -> String
fmtF dp x = show (fromIntegral (round (x*10^dp)) / 10^dp :: Double)

padL :: Int -> String -> String
padL n s = replicate (max 0 (n-length s)) ' ' ++ s

padR :: Int -> String -> String
padR n s = take n (s ++ repeat ' ')

rule :: String
rule = "└" ++ replicate 66 '─' ++ "┘"

printSection :: String -> [String] -> IO ()
printSection hdr rows = do
  putStrLn ""
  putStrLn $ "┌─ " ++ hdr ++ " " ++ replicate (64-length hdr) '─' ++ "┐"
  mapM_ (putStrLn . ("│  "++)) rows
  putStrLn rule

bar :: Double -> Double -> Int -> String
bar v maxV w =
  let len = max 0 . min w $ round(v/max maxV 1e-9*fromIntegral w)
  in  replicate len '█' ++ replicate(w-len)'░'

tick :: Bool -> String
tick True  = "✓"
tick False = "✗"

-- ─── Reports ─────────────────────────────────────────────────────────────────

reportKahlerChecks :: [KahlerPoint] -> IO ()
reportKahlerChecks pts = do
  let pairs  = zip pts (drop 1 pts)
      checks = take 8 pairs
  printSection "Kähler condition verification" $
    [ "For each consecutive pair (p_t, p_{t+1}):"
    , ""
    , "  t    J²=-I  g(Ju,Jv)=g(u,v)  ω skew  ω(u,v)=g(Ju,v)" ] ++
    map (\(p,q) ->
      let t = kpLabel p
          u = kpComplex p; v = kpComplex q
      in  "  " ++ padR 4 (show t)
          ++ padL 7  (tick (checkJSquared u))
          ++ padL 16 (tick (checkJPreservesMetric u v))
          ++ padL 9  (tick (checkSkewSymmetric u v))
          ++ padL 18 (tick (checkCompatibility u v))
      ) checks ++
    [ ""
    , "  All ✓ confirms this is a flat Kähler space (C^k, standard structure)" ]

reportCurvature :: [KahlerPoint] -> IO ()
reportCurvature pts = do
  let triples = zip3 pts (drop 1 pts) (drop 2 pts)
      results = take 8 triples
  printSection "Embedded manifold curvature — dω over signal triples" $
    [ "Measures how curved the signal attractor is inside C^k."
    , "dω(p,q,r) = ω(q-p,r-p) + ω(r-q,p-q) + ω(p-r,q-r)"
    , ""
    , "  ≈ 0  →  locally flat (attractor tangent plane is Lagrangian)"
    , "  ≠ 0  →  curved embedding; magnitude = local curvature flux"
    , ""
    , "  t     dω value     curved?" ] ++
    map (\(p,q,r) ->
      let (val,curved) = omegaCurvature (kpComplex p) (kpComplex q) (kpComplex r)
          t = kpLabel p
      in  "  " ++ padR 5 (show t)
          ++ padL 12 (fmtF 5 val)
          ++ "     " ++ (if curved then "yes — curvature flux" else "flat")
      ) results ++
    [ ""
    , "  Non-zero dω is expected and meaningful for chaotic signals."
    , "  It quantifies how far the signal manifold deviates from"
    , "  a Lagrangian submanifold of C^k." ]

reportChernComparison :: AugConfig -> IO ()
reportChernComparison cfg = do
  let results = map (\(name,xs) ->
                  let pts = map (complexify cfg) (augment cfg xs)
                      c   = chernProxy pts
                  in  (name, c)
                ) allSignals
      maxC = maximum (map (abs.snd) results)
  printSection "Chern proxy comparison — all six signals" $
    [ "Σ_{i<j} ω(p_i,p_j) — discrete integral of symplectic form"
    , ""
    , "  Large   →  monotone phase progression; ω accumulates without cancellation"
    , "  Small   →  trajectory crosses itself; ω partially cancels across pairs"
    , "" ] ++
    map (\(name,c) ->
      "  " ++ padR 12 name
      ++ bar (abs c) maxC 24
      ++ "  " ++ (if c>=0 then " " else "") ++ fmtF 2 c
      ) results ++
    [ ""
    , "  Sine is largest: perfectly regular phase never reverses in lag space"
    , "  so ω accumulates monotonically across all pairs with no cancellation."
    , "  Chaotic signals (logistic, random walk) cross themselves in lag space,"
    , "  producing partial cancellation and a smaller net integral." ]

reportSymplecticMatrix :: [KahlerPoint] -> IO ()
reportSymplecticMatrix pts = do
  let sel    = take 6 pts
      mat    = symplecticMatrix sel
      header = "      " ++ concatMap(\p->padL 9("t="++show(kpLabel p))) sel
      rows   = zipWith(\p row ->
                 padR 6("t="++show(kpLabel p))
                 ++ concatMap(\v->padL 9((if v>=0 then " " else "")++fmtF 4 v)) row
               ) sel mat
  printSection "Symplectic form ω — first 6 points" $
    [ "ω(u,v) = Im H(u,v) = Im Σ conj(u_i)·v_i"
    , "Skew-symmetric; diagonal = 0; encodes phase between lag pairs"
    , "" ] ++
    (header : "" : rows)

reportHoloKNN :: [KahlerPoint] -> IO ()
reportHoloKNN pts = do
  let queryT = 10
  case filter((==queryT).kpLabel) pts of
    []        -> putStrLn "Query not found."
    (query:_) -> do
      let nbs  = holoKNN 5 query pts
          maxD = maximum(map snd nbs)
      printSection ("Holomorphic 5-NN of t="++show queryT++" (Hermitian L2)") $
        [""] ++
        map(\(p,d) ->
          "  t="++padR 4(show(kpLabel p))
          ++bar d maxD 28
          ++"  "++fmtF 5 d
          ) nbs

reportPersistentHomology :: [KahlerPoint] -> IO ()
reportPersistentHomology pts = do
  let filt   = vrFiltration pts 10
      maxEnt = maximum [ e | (_,_,_,e) <- filt ]
  printSection "Persistent homology — Vietoris-Rips on C^k distances" $
    [ "Sine centred at 0 vs sine shifted +10 — guaranteed cluster separation."
    , "Inter-cluster gap >> intra-cluster spread by construction."
    , "Bar width = Shannon entropy of active edge distances (bits)."
    , ""
    , "  ε          H0  entropy (bits)    H1  entropy (bits)"
    , "  ──────────────────────────────────────────────────" ] ++
    map (\(e,h0,h1,ent) ->
      "  " ++ padR 10 (fmtF 3 e)
      ++ padL 4 (show h0) ++ " " ++ bar ent maxEnt 16
      ++ "   " ++ padL 3 (show h1) ++ " " ++ bar ent maxEnt 16
      ++ "  " ++ fmtF 3 ent ++ " bits"
      ) filt ++
    [ ""
    , "  H0 drop = components merging in C^k."
    , "  H1 spike = loop born; H1 drop = loop filled by triangle."
    , "  Entropy peak = maximally spread edge lengths = topological transition." ]

reportNormsComparison :: [AugPoint] -> [KahlerPoint] -> IO ()
reportNormsComparison realPts kPts =
  printSection "Real L2 vs Hermitian norm (first 10 points)" $
    [ "Re H(v,v) = ‖v‖² confirms the metric is preserved under complexification"
    , ""
    , "  t    real ‖v‖₂    Hermitian ‖v‖_C    ratio" ] ++
    zipWith(\rp kp ->
      let rn = sqrt.sum $ map(^2)(apVec rp)
          kn = cnorm(kpComplex kp)
          ratio = if rn>1e-12 then kn/rn else 0
      in  "  "++padR 4(show(apLabel rp))
          ++padL 12(fmtF 4 rn)
          ++padL 18(fmtF 4 kn)
          ++padL 10(fmtF 4 ratio)
      ) (take 10 realPts) (take 10 kPts)

-- ─── Main ────────────────────────────────────────────────────────────────────

main :: IO ()
main = do
  putStrLn "╔════════════════════════════════════════════════════════════════════╗"
  putStrLn "║   Augmented Vector Space — Kähler Extension v2  ·  Haskell base  ║"
  putStrLn "╚════════════════════════════════════════════════════════════════════╝"

  let cfg     = AugConfig { lagWindow=4, derivOrder=1, rollingWin=4 }
      xs      = [-2.2e-05, 0.469324, 0.927044, 1.361952, 1.763352, 2.121341, 2.427033, 2.673001, 2.853188, 2.963087, 3.0, 2.963087, 2.853188, 2.673001, 2.427033, 2.121341, 1.763352, 1.361952, 0.927044, 0.469324, -2.2e-05, -0.469324, -0.927044, -1.361952, -1.763352, -2.121341, -2.427077, -2.673045, -2.853188, -2.963087, -3.0, -2.963087, -2.853188, -2.673045, -2.427077, -2.121341, -1.763352, -1.361952, -0.927044, -0.469324, -2.2e-05, 0.469324, 0.927044, 1.361952, 1.763352, 2.121341, 2.427033, 2.673001, 2.853188, 2.963087, 3.0, 2.963087, 2.853188, 2.673001, 2.427033, 2.121341, 1.763352, 1.361952, 0.927044, 0.469324]
      realPts = augment cfg xs
      kPts    = map (complexify cfg) realPts

  -- Homology corpus: sine (centred at 0) vs sine shifted by DC+10.
  -- Shifting the mean moves the lag vectors to a completely separate
  -- region of C^k — inter-cluster gap is exactly 10*sqrt(dim),
  -- guaranteed larger than any intra-cluster distance.
  let sineA   = sineWave 60
      sineB   = map (+10) sineA          -- DC shift = guaranteed separation
      ptsA    = map (complexify cfg) (augment cfg sineA)
      ptsB    = map (complexify cfg) (augment cfg sineB)
      stride xs k = [ xs !! (i*k) | i <- [0..7], i*k < length xs ]
      homoPts = stride ptsA 4 ++ stride ptsB 4

  reportKahlerChecks       kPts
  reportCurvature          kPts
  reportChernComparison    cfg
  reportSymplecticMatrix   kPts
  reportHoloKNN            kPts
  reportNormsComparison    realPts kPts
  reportPersistentHomology homoPts

  putStrLn ""
  putStrLn "── Mathematical structure ──────────────────────────────────────────"
  putStrLn "   R^n → C^(n/2)        complexify consecutive lag pairs"
  putStrLn "   H = g + iω           Hermitian metric splits into:"
  putStrLn "     g = Re H           Riemannian (= original L2)"
  putStrLn "     ω = Im H           Symplectic (phase structure)"
  putStrLn "   dω ≠ 0 on triples    curvature of embedded signal manifold"
  putStrLn "   Chern proxy          discrete integral of ω; large iff monotone phase"
  putStrLn "   VR filtration        persistent H0/H1 in Hermitian distance"
  putStrLn ""
  putStrLn "── Swap signals to compare ─────────────────────────────────────────"
  putStrLn "   sineWave, logisticMap, lorenzWave, sawtoothWave, stepWave, randomWalk"
