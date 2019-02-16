{-# LANGUAGE ForeignFunctionInterface #-}
module CodeBlock (plugin) where

-- This plugin reformats C code after processing by pandoc.

import Foreign
import Foreign.C
import System.IO.Unsafe as Unsafe
import Control.Monad.State (get)
import Network.Gitit.Interface

plugin :: Plugin
plugin = mkPageTransformM fixBlock

fixBlock :: Block -> PluginM Block
fixBlock (CodeBlock attr@(_,classes,_) s) | "literatec" `elem` classes = do
  ctx <- get
  return $ CodeBlock attr (codeSubst s (pgPageName $ ctxLayout ctx))
fixBlock x = return x

-- the actual transformation is done in C, Haskell is just too hard...
foreign import ccall unsafe "codeblock"
     codeblock_c :: CString -> CString -> IO CString

-- from http://rosettacode.org/wiki/Call_a_foreign-language_function#Haskell
codeSubst :: String -> String -> String
codeSubst s pname = Unsafe.unsafePerformIO $
      withCString s $ \cs -> do
      withCString pname $ \cname -> do
	s2 <- codeblock_c cs cname
        s2_hs <- peekCString s2 -- marshall the C string called s2 into a Haskell string named s2_hs
	free s2
	return s2_hs -- cs is automatically freed by withCString once done
