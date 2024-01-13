use proc_macro::{Span, TokenStream};
use syn::{
    parse_macro_input, parse_quote, Error, GenericParam, Generics, Item,
    TypeParamBound,
};

/// Adds a trait bound to every generic parameter in an impl block or struct
/// definition
fn add_trait_bound(
    mut item: Item,
    traits: &[TypeParamBound],
) -> Result<Item, Error> {
    if let Item::Impl(ref mut impl_item) = item {
        traits
            .iter()
            .for_each(|t| append_generic(&mut impl_item.generics, t));
        Ok(item)
    } else if let Item::Struct(ref mut struct_item) = item {
        traits
            .iter()
            .for_each(|t| append_generic(&mut struct_item.generics, t));
        Ok(item)
    } else {
        Err(Error::new(
            Span::call_site().into(),
            "macro must be called on impl block or struct definition",
        ))
    }
}

/// Add a bound to every type parameter
fn append_generic(generics: &mut Generics, picked_trait: &TypeParamBound) {
    for param in &mut generics.params {
        if let GenericParam::Type(ref mut type_param) = *param {
            type_param.bounds.push(picked_trait.clone());
        }
    }
}

fn feature_check() -> Result<(), Error> {
    if cfg!(all(feature = "sync", feature = "tokio"))
        || cfg!(not(any(feature = "sync", feature = "tokio")))
    {
        Err(Error::new(
            Span::call_site().into(),
            "One and only one sync/async feature must be enabled",
        ))
    } else {
        Ok(())
    }
}

/// Adds the trait bounds required for either sync or async writing
///
/// These bounds are applied to all generic parameters in either an impl block or
/// struct definition.
///
/// Traits are chosen according to which features are enabled.
#[proc_macro_attribute]
pub fn write_bound(_: TokenStream, item: TokenStream) -> TokenStream {
    if let Err(e) = feature_check() {
        return Error::into_compile_error(e).into();
    }
    let write_trait = if cfg!(feature = "sync") {
        vec![parse_quote!(::std::io::Write)]
    } else if cfg!(feature = "tokio") {
        vec![
            parse_quote!(::tokio::io::AsyncWriteExt),
            parse_quote!(::std::marker::Unpin),
        ]
    } else {
        unreachable!()
    };
    let item = add_trait_bound(parse_macro_input!(item as Item), &write_trait);
    match item {
        Ok(item) => quote::quote!(#item).into(),
        Err(e) => Error::into_compile_error(e).into(),
    }
}
