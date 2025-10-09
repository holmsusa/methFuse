# bino_div_R
test_that("bino_div works as expected", {
  expect_equal(.Call("bino_div_R", c(4, 4, 4, 0, -3), c(2, 0, -1, 5, 1)),
               c(2, 0.5, 0.5, 0, -3))
  expect_error(.Call("bino_div_R", 1:3, 1:2))
  expect_error(.Call("bino_div_R", "a", "b"))
  expect_equal(.Call("bino_div_R", numeric(0), numeric(0)), numeric(0))
})


# bino_xlogy_R
test_that("bino_xlogy works as expected", {
  expect_equal(
    .Call("bino_xlogy_R", c(2, 2, 2, 2, -1),
          c(1, exp(1), 0, -1, exp(2))),
    c(0, 2, 0, 0, -2)
  )
  expect_error(.Call("bino_xlogy_R", 1:3, 1:2))
  expect_error(.Call("bino_xlogy_R", "a", "b"))
  expect_equal(.Call("bino_xlogy_R", numeric(0), numeric(0)), numeric(0))
})
