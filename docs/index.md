# Documentation

Start with [getting started](quickstart.md) for installation and a first geometric calculation. Use the guides below to choose an object, interpret its result and check its mathematical scope. The [API reference](api.md) supplies signatures and result types after that choice.

## Objects, calculations and results

| Question | Read |
| --- | --- |
| Which geometric object should represent my equations? | [Conventions](conventions.md), [geometry](geometry.md), [schemes](schemes.md) |
| How do I read dimensions, multiplicities and certificates? | [Capabilities and limits](capabilities.md), [result interpretation](results.md), [mathematical methods](methods.md) |
| How do I compute an example from start to finish? | [Worked applications](worked_applications.md), [application calculations](applications.md) |
| Which coefficient field and backend should I use? | [Backend installation](backends.md), [field support](coefficient_fields_support.md), [coefficient fields](coefficient_fields.md), [number fields](number_fields.md), [finite fields](finite_fields.md) |

## Algebraic geometry

Choose [components](components.md) for radicals and primary structure, [local geometry](local_geometry.md) for pointwise invariants, and [normalization](normalization.md) or [curves](curves.md) for birational models. [Modules](modules.md) and [resolutions](resolutions.md) describe relations and homological invariants.

For maps and images, begin with [geometric constructions](geometric_constructions.md), then use [constructive Chevalley projection](constructive_chevalley.md) for exact locally closed images. [Absolute factorization](absolute_factorization.md) and [constant fields and geometric primality](function_field_constants.md) explain how changing the ground field affects components.

## Parameters and specialization

A generic calculation over a function field can omit exceptional fibers. Use [nonlinear parameter families](nonlinear_families.md) and [coupled fiber stratification](comprehensive_fibers.md) when those fibers matter.

[Gröbner covers](groebner_covers.md) introduce polynomial bases on parameter strata. Continue with [prime parameter ideals](parameter_ideals.md), [quotient localizations and verified Gröbner construction](quotient_localization_and_groebner_backend.md), and [quadratic constants and localized linear covers](quadratic_fields_and_linear_covers.md) for the supported arithmetic and certificate contracts.

## Function decomposition and real domains

The [decomposition guide](guide.md) distinguishes rational functions, implicit relations and selected real branches. It explains search completeness, component reconstruction and the optional semialg adapter. [Package ownership](package_ownership.md) shows which computations belong to alggeom, algroots and semialg.

## Validation and development

Use [input limitations](known_issues.md) to recognize unsupported cases, [testing](testing.md) for independent mathematical expectations and installed examples, and [certificate audits](certificate_audits.md) for replay obligations. [Benchmarks](benchmarks.md) describe reproducible performance measurements; [release checks](releasing.md) collect the steps needed before publication.

The [mathematical references](references.md) and [bibliography](../REFERENCES.md) give the foundations behind the implementations. The [API reference](api.md) links each public object to its defining source.
