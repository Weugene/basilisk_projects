{-
Copyright (C) 2012 John Lenz <lenz@math.uic.edu>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.  Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.  Neither the name of John Lenz nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-}

module BibtexGitit where

-- For documentation, see
-- * http://blog.wuzzeb.org/posts/2012-06-19-bibtex-and-pandoc-2.html
-- * http://blog.wuzzeb.org/posts/2012-06-26-bibtex-and-gitit.html

import Control.Applicative ((<$>))
import Control.Exception (throw)
import Control.Monad (liftM)
import Data.Maybe (fromMaybe)
import Data.Char (toLower, isAlpha)
import qualified Data.List as L
import qualified Data.List.Split as S
import System.FilePath (replaceExtension)
import System.Environment (getArgs)
import qualified Text.ParserCombinators.Parsec as P
import Text.ParserCombinators.Parsec ((<|>))
import Text.Pandoc

import Network.Gitit.Interface

data Bibtex = Bibtex String [(String,String)]

bibParser :: P.Parser [Bibtex]
bibParser = do
  x <- P.sepEndBy bibEntry P.spaces
  P.eof
  return x

bibEntry :: P.Parser Bibtex
bibEntry = do
  P.char '@'
  P.many $ P.noneOf "{"
  P.char '{'
  name <- P.many $ P.noneOf ","
  P.many $ P.oneOf " \t\n,"
  attrs <- P.sepEndBy bibAttr $ P.many1 $ P.oneOf " \t\n,"
  P.char '}'
  return $ Bibtex name attrs

bibAttr :: P.Parser (String, String)
bibAttr = do
  key <- P.many (P.letter <|> P.digit)
  P.spaces
  P.char '='
  P.spaces
  P.char '{'
  val <- bibVal
  P.char '}'
  return (map toLower key, val)

bibVal :: P.Parser String
bibVal = liftM concat $ P.many1 (bibValMatched <|> (liftM (:[]) (P.noneOf "{}")))

bibValMatched :: P.Parser String
bibValMatched = P.between (P.char '{') (P.char '}') bibVal

renderEntries :: [Bibtex] -> Block
renderEntries lst = DefinitionList $ map display lst'
    where lst' = L.sortBy (\(Bibtex a _) (Bibtex b _) -> compare a b) lst
          display (Bibtex key b) = ([Strong [Str key]], [[Plain $ renderEntry key b]])

type BibtexAttr = [(String, String)]

render1 :: BibtexAttr -> String -> Inline
render1 b s = case lookup s b of
    Just x  -> Str x
    Nothing -> Str ""

articleLink :: String -> BibtexAttr -> Inline
articleLink s b = case lookup s b of
                  Just x  -> Link [Str "article"] (x, [])
                  Nothing -> Str ""

mrNumber :: BibtexAttr -> Inline
mrNumber b = case lookup "mrnumber" b of
                Just x -> mkURL x
                Nothing -> Str ""
  where mkURL x | length x > 2 = Link [Str "MathSciNet"] (mathSciNet ++ mrNum x, [])
        mkURL _ = Str ""
        mrNum = dropWhile isAlpha . head . words
        mathSciNet = "http://www.ams.org/mathscinet/search/publdoc.html?pg1=MR&s1="

arxiv :: BibtexAttr -> Inline
arxiv b = case lookup "arxiv" b of
                 Just x  -> mkURL x
                 Nothing -> Str ""
 where mkURL x = Link [Str "arXiv"] (url ++ dropWhile isAlpha x, [])
       url = "http://arxiv.org/abs/"

expandTex :: String -> String
expandTex ('\\':a:'{':b:'}':xs) = expandTex ('\\':a:b:xs)
expandTex ('\\':'\'':a:xs) = a' : expandTex xs
   where a' = case a of
               'a' -> 'á'
               'e' -> 'é'
               'o' -> 'ó'
               _   -> a
expandTex ('\\':'H':'o':xs) = 'ő' : expandTex xs
expandTex ('\\':'"':a:xs) = a' : expandTex xs
   where a' = case a of
               'a' -> 'ä'
               'e' -> 'ë'
               'o' -> 'ö'
               _   -> a
expandTex (a:xs) = a : expandTex xs
expandTex [] = []

prettyAuthor :: String -> String
prettyAuthor x = L.intercalate ", " $ map fixOne $ S.splitOn " and" x
  where fixOne s = case S.splitOn "," s of
                    []     -> ""
                    [a]    -> a
                    (f:xs) -> concat xs ++ " " ++ f

renderEntry :: String -> BibtexAttr -> [Inline]
renderEntry name b = raw ++ entries
  where raw = [(RawInline "html" $ "<a name=\"" ++ name ++ "\"></a>")]

        entries = L.intersperse (Str ", ") $ filter (not . isEmptyStr)
            [ mapInline (prettyAuthor . expandTex) $ render1 b "author"
            , mapInline (\a -> "\"" ++ a ++ "\"") $ render1 b "title"
            , Emph [render1 b "journal"]
            , render1 b "year"
            , mrNumber b
            , articleLink "url" b
            , articleLink "url2" b
            , arxiv b
            ]

        mapInline f (Str s) = Str $ f s
        mapInline _ x = x
        isEmptyStr (Str "") = True
        isEmptyStr _        = False

transformBlock :: Block -> Block
transformBlock (CodeBlock (_, classes, namevals) contents) | "bib" `elem` classes =
        case P.parse bibParser "" contents of
           Left err -> BlockQuote [Para [Str $ "Error parsing bib data: " ++ show err]]
           Right x -> renderEntries x
transformBlock x = x

plugin :: Plugin
plugin = mkPageTransform transformBlock
